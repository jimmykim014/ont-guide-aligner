"""Core counting logic for ont-guide-aligner.

Each Oxford Nanopore read counts toward at most ONE guide. How a read is assigned:

  1. The vector sequence immediately upstream and downstream of the guide (``flank_len`` bp
     each, i.e. flank | guide | flank) is learned per sample from the first reads: the most
     common flanks around exact guide matches, scanning until ``flank_hits`` exact matches or
     ``sample_size`` reads.
  2. Both orientations of every read are searched for those flanks anywhere in the read. Each
     flank occurrence marks where a guide should start, so guides are found whether the read
     is full length, truncated at either end, or longer than expected (ONT indels shift the
     guide's position from read to read).
  3. The guide-length window at that spot must match a library guide with at most 1 mismatch
     (default) or exactly (``exact_match``, e.g. base-editor libraries). A 1-mismatch sequence
     that is 1 base from two different guides is never used; a read whose only such match is
     one of those is dropped as ambiguous.
  4. Matches are ranked by (intact flanks, exact?): both flanks beat one, exact beats 1
     mismatch. The read goes to the single best guide; a tie between different guides (mostly
     two amplicons joined in one read) is dropped as ambiguous.

Both orientations are handled in one pass, so a read can never be counted twice.

Indel rescue (``allow_indel``, on by default): a read with no match in steps 3-4 gets a second
chance where both flanks are found glen +/- 1 bp apart. The bases between them may be one
insertion or deletion from a library guide; the read is counted if exactly one guide fits,
else dropped as ambiguous.

Flank learning (step 1) runs here in Python on the first reads; steps 2-4 and the indel rescue
run in the compiled Rust core (``ont_guide_aligner._core``, src/lib.rs), which is always built
with the package.
"""

import csv
import gzip
import math
import os
import signal
import zlib
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool

try:
    from ._core import count_file as _count_file
except ImportError as e:  # only possible when running from an unbuilt source checkout
    raise ImportError(
        "ont-guide-aligner's compiled core (ont_guide_aligner._core) is missing. Reinstall with "
        "pip, or build a source checkout with: maturin develop --release") from e

__all__ = [
    "LibraryError",
    "read_library",
    "count_fastqs",
    "qc_stats",
    "sample_name",
]

_RC = str.maketrans("ACGTN", "TGCAN")
_NO_MAX_LEN = 2**32 - 1  # "no maximum read length"; fits the Rust core's usize on every target


class LibraryError(ValueError):
    """The guide library can't be used as given (missing columns or values, non-ACGT bases,
    duplicate IDs or sequences, mixed lengths)."""


def revcomp(s):
    return s.translate(_RC)[::-1]


# Header names recognised for each library column, compared lower-case with spaces,
# underscores and other punctuation removed ("sgRNA Sequence" -> "sgrnasequence").
# Names are matched whole, so GeCKO's "gene_id" is a gene column, not an ID ("UID" is the ID).
_ID_NAMES = {"guideid", "guide", "guidename", "id", "name", "sgrna", "sgrnaid", "sgrnaname", "uid"}
_SEQ_NAMES = {"sequence", "seq", "guidesequence", "sgrnasequence", "spacer", "spacersequence",
              "sgrnatargetsequence", "targetsequence", "sgrnaseq", "protospacer",
              "protospacersequence"}
_GENE_NAMES = {"gene", "geneid", "genename", "genesymbol", "symbol", "target", "targetgene",
               "targetgenesymbol"}


def _norm_header(name):
    return "".join(c for c in name.lower() if c.isalnum())


def _is_dna(s):
    return bool(s) and set(s.upper()) <= set("ACGT")


def read_library(path):
    """Read a guide library: tab- or comma-separated, one guide per row with an ID, a sequence
    and a gene. Returns (ids, seqs, genes, guide_length).

    The first row is taken as a header unless its second field is already a DNA sequence (or
    looks like a mistyped one: no recognised column name, and letters only of the same length
    as the next row's sequence). Columns are found by header name (e.g. guide_id / guide,
    sequence / sgRNA Sequence, gene; see _ID_NAMES etc.); when the sequence column's name
    isn't recognised, the columns are taken in order: ID, sequence, gene.

    Raises LibraryError for a header without a gene column, a row without an ID, sequence or
    gene, a sequence with characters other than A/C/G/T, duplicated guide IDs or sequences
    (their counts would be merged or unassignable), or guides of different lengths."""
    with open(path, "rb") as f:
        bom = f.read(2)
    # Excel's "Unicode Text" export is UTF-16 with a BOM; otherwise UTF-8 (-sig: drop a BOM)
    encoding = "utf-16" if bom in (b"\xff\xfe", b"\xfe\xff") else "utf-8-sig"
    with open(path, newline="", encoding=encoding) as f:
        sample = f.read(4096)
        f.seek(0)
        delim = "\t" if "\t" in sample else ","
        reader = csv.reader(f, delimiter=delim)
        rows = [(reader.line_num, r) for r in reader if any(x.strip() for x in r)]
    if not rows:
        raise LibraryError(f"{path}: the file is empty.")

    first = [x.strip() for x in rows[0][1]]
    has_header = not (len(first) > 1 and _is_dna(first[1]))
    if has_header and len(first) > 1 and len(rows) > 1:
        known = _ID_NAMES | _SEQ_NAMES | _GENE_NAMES
        nxt = [x.strip() for x in rows[1][1]]
        # A first guide with a bad sequence, not a header: let the row checks below report it
        if (not any(_norm_header(x) in known for x in first) and first[1].isalpha()
                and len(nxt) > 1 and len(first[1]) == len(nxt[1])):
            has_header = False
    id_col, seq_col, gene_col = 0, 1, 2
    by_name = False
    if has_header:
        names = [_norm_header(x) for x in first]
        find = lambda wanted: next((i for i, n in enumerate(names) if n in wanted), None)
        if find(_SEQ_NAMES) is not None:
            by_name = True
            seq_col = find(_SEQ_NAMES)
            id_col = find(_ID_NAMES)
            gene_col = find(_GENE_NAMES)
            if id_col is None or id_col == seq_col:
                raise LibraryError(
                    f"{path}: found the sequence column ({first[seq_col]!r}) but no guide ID "
                    f"column in the header {first}. Name it guide_id.")
            if gene_col is None or gene_col in (id_col, seq_col):
                raise LibraryError(
                    f"{path}: found the guide ID and sequence columns but no gene column in the "
                    f"header {first}. Add one named gene.")
        rows = rows[1:]
    if not rows:
        raise LibraryError(f"{path}: no guides found after the header row.")

    ids, seqs, genes = [], [], []
    for line, r in rows:
        r = [x.strip() for x in r]
        gid = r[id_col] if id_col < len(r) else ""
        seq = r[seq_col].upper() if seq_col < len(r) else ""
        gene = r[gene_col] if gene_col < len(r) else ""
        if not gid or not seq or not gene:
            what = "guide ID" if not gid else "guide sequence" if not seq else "gene"
            raise LibraryError(
                f"{path}, line {line}: missing the {what} (row: {delim.join(r)!r}). Each row "
                "needs a guide ID, a sequence and a gene.")
        if not _is_dna(seq):
            bad = sorted(set(seq) - set("ACGT"))
            raise LibraryError(
                f"{path}, line {line}: guide {gid!r} has a sequence with characters other than "
                f"A, C, G, T ({', '.join(bad)}): {seq!r}."
                + ("" if by_name else " Columns were taken in order (ID, sequence, gene); if "
                   "yours are in a different order, name them in a header row: guide_id, "
                   "sequence, gene."))
        ids.append(gid)
        seqs.append(seq)
        genes.append(gene)

    dup_ids = [k for k, v in Counter(ids).items() if v > 1]
    if dup_ids:
        raise LibraryError(f"{len(dup_ids)} guide IDs are used more than once (e.g. {dup_ids[0]!r}). "
                           "Give every guide a unique ID first.")
    dup_seqs = [k for k, v in Counter(seqs).items() if v > 1]
    if dup_seqs:
        raise LibraryError(f"{len(dup_seqs)} sequences appear under more than one guide "
                           f"(e.g. {dup_seqs[0]}). Remove or merge them first.")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise LibraryError(f"Guides must all be the same length; found lengths {sorted(lengths)}.")
    return ids, seqs, genes, lengths.pop()


def _seq_lines(f, path):
    """Yield the sequence line of each 4-line FASTQ record, checking the record structure the
    same way the Rust core does (blank lines between records are allowed)."""
    pos = 0
    seq_len = 0
    for line in f:
        if pos == 0:
            if not line.strip():
                continue
            if not line.startswith("@"):
                raise RuntimeError(f"{path}: not a valid FASTQ file: a read's first line doesn't "
                                   "start with '@' (it may be truncated or corrupted)")
        elif pos == 1:
            seq = line.strip().upper()
            seq_len = len(seq)
            yield seq
        elif pos == 2 and not line.startswith("+"):
            raise RuntimeError(f"{path}: not a valid FASTQ file: the line after a read's "
                               "sequence doesn't start with '+' (it may be truncated or corrupted)")
        elif pos == 3 and len(line.strip()) != seq_len:
            raise RuntimeError(f"{path}: not a valid FASTQ file: a read's quality line is not the "
                               "same length as its sequence (it may be truncated or corrupted)")
        pos = (pos + 1) % 4


def iter_reads(path):
    """Yield read sequences from a FASTQ or FASTQ.gz file. Only flank learning reads files in
    Python, and it stops after the first reads, so the gzip module is fast enough."""
    with open(path, "rb") as f:
        is_gz = f.read(2) == b"\x1f\x8b"  # by content, so .FASTQ.GZ, .bgz etc. work
    opener = gzip.open if is_gz else open
    try:
        # only sequence lines are used, so odd bytes in headers mustn't stop the run
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            yield from _seq_lines(f, path)
    except (EOFError, zlib.error, gzip.BadGzipFile) as e:
        raise RuntimeError(f"{path}: can't decompress ({e}); the file may be truncated or "
                           "corrupted") from e


def detect_flanks(path, index, glen, flen, sample_size, max_hits, min_len=0, max_len=None):
    """Learn the vector flanks from exact guide matches in the first reads.

    Scans up to `sample_size` reads, stopping early once `max_hits` exact guide matches have
    been seen. Reads outside min_len..max_len are skipped, as they are when counting. Returns a
    dict: up, down (most common flank on each side, in guide orientation; None if a side was
    never observed), up_frac, down_frac (their shares) and hits."""
    up, down = Counter(), Counter()
    hits = 0
    for n, read in enumerate(iter_reads(path)):
        if n >= sample_size or hits >= max_hits:
            break
        L = len(read)
        if L < min_len or (max_len is not None and L > max_len):
            continue
        for seq in (read, revcomp(read)):
            for o in range(L - glen + 1):
                if seq[o:o + glen] in index:
                    hits += 1
                    if o >= flen:
                        up[seq[o - flen:o]] += 1
                    if o + glen + flen <= L:
                        down[seq[o + glen:o + glen + flen]] += 1
    flanks = {"up": None, "down": None, "up_frac": 0.0, "down_frac": 0.0, "hits": hits}
    if up and down:
        (u, uc), (d, dc) = up.most_common(1)[0], down.most_common(1)[0]
        flanks.update(up=u, down=d, up_frac=uc / sum(up.values()), down_frac=dc / sum(down.values()))
    return flanks


def gini(values):
    v = sorted(values)
    n, total = len(v), sum(v)
    if n == 0 or total == 0:
        return float("nan")
    return sum((2 * (i + 1) - n - 1) * x for i, x in enumerate(v)) / (n * total)


def percentile(values, q):
    """q-th percentile (0-100) with linear interpolation between ranks, as numpy.percentile."""
    v = sorted(values)
    if not v:
        return float("nan")
    pos = (len(v) - 1) * q / 100
    lo = int(pos)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


def skew_ratio(counts):
    """90th / 10th percentile guide count, or None when the 10th percentile is 0 (many
    zero-count guides, already shown by Zero-Count Guides and Guide Representation)."""
    if not counts:
        return None
    p10 = percentile(counts, 10)
    return percentile(counts, 90) / p10 if p10 else None


def qc_stats(counts, stats):
    """Library-representation QC for one sample, from its guide counts and the stats dict
    count_fastqs returned for it. These are the columns of <output>.QC.stats.txt; the Gini
    coefficient is "NA" when no reads were counted, the skew ratio when the 10th-percentile
    count is 0."""
    n = len(counts)
    total_counts = sum(counts)
    zero = sum(1 for c in counts if c == 0)
    skew = skew_ratio(counts)
    g = gini(counts)
    return {
        "Total Guides": n,
        "Zero-Count Guides": zero,
        "Guide Representation (%)": round(100 * (n - zero) / n, 2) if n else 0.0,
        "Total Counts": total_counts,
        "Average Count": round(total_counts / n, 2) if n else 0.0,
        "Max Count": max(counts, default=0),
        "Gini Coefficient": "NA" if math.isnan(g) else round(g, 3),  # NA: no reads counted
        "Skew Ratio (P90/P10)": round(skew, 2) if skew is not None else "NA",
        "Total Reads": stats["total_reads"],
        "Reads Assigned": stats["assigned_reads"],
        "Percent Assigned (%)": stats["pct_assigned"],
    }


def count_sample(job):
    (path, seqs, glen, flen, sample_size, flank_hits, min_len, max_len, exact_only, threads,
     allow_indel) = job
    index = {s: i for i, s in enumerate(seqs)}
    fl = detect_flanks(path, index, glen, flen, sample_size, flank_hits, min_len, max_len)
    if not fl["hits"]:
        raise RuntimeError(f"{path}: no exact guide matches in the first {sample_size} reads; "
                           "check the library and that this FASTQ belongs to it.")
    up, down = fl["up"], fl["down"]
    if up is None:
        raise RuntimeError(
            f"{path}: can't learn vector flanks on both sides of the guide. ont-guide-aligner "
            "needs reads that span the guide plus some vector on each side (ONT amplicon reads); "
            "reads that end at the guide, such as short Illumina reads, aren't supported.")

    counts, t = _count_file(path, seqs, glen, up, down, min_len, max_len, exact_only, threads,
                            allow_indel)
    total, assigned, both_flank = t["total"], t["assigned"], t["both_flank"]

    stats = {
        "total_reads": total,
        "assigned_reads": assigned,
        "pct_assigned": round(100 * assigned / total, 2) if total else 0.0,
        "assigned_both_flanks": both_flank,
        "assigned_one_flank": assigned - both_flank,
        "assigned_1_mismatch": t["mismatch"],
        "assigned_indel": t["indel"] if allow_indel and not exact_only else "NA",
        "ambiguous_dropped": t["ambiguous"],
        # reads skipped by an optional --min/--max-read-len limit count as unmatched
        "unmatched": t["unmatched"] + t["length_filtered"],
        "upstream_flank": f"{up} ({100 * fl['up_frac']:.0f}%)",
        "downstream_flank": f"{down} ({100 * fl['down_frac']:.0f}%)",
        # for the CLI's flank warnings; keys starting with "_" aren't written to the stats files
        "_flank_hits": fl["hits"],
        "_up_frac": fl["up_frac"],
        "_down_frac": fl["down_frac"],
        "_flanks": (up, down),
    }
    return list(counts), stats


def sample_name(path):
    """Sample name from a FASTQ path: the file name without .gz/.bgz/.fastq/.fq (any case)."""
    name = os.path.basename(path)
    for ext in (".gz", ".bgz", ".fastq", ".fq"):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
    return name


def count_fastqs(paths, seqs, *, exact_match=False, flank_len=6, sample_size=100000, flank_hits=2000,
                 min_len=0, max_len=None, jobs=None, threads=None, allow_indel=True,
                 skip_failed=False, on_done=None):
    """Count each FASTQ in `paths` against guide sequences `seqs` (all the same length, unique).

    The caller owns guide IDs and genes; results are lists aligned with `seqs`. Returns one
    (counts, stats) tuple per path, in input order, where stats is the dict written to
    <output>.alignment.stats.txt (pass both to qc_stats for <output>.QC.stats.txt). Samples run
    in parallel processes (`jobs`); `threads` splits each multi-member .gz into pieces counted
    in parallel. `allow_indel` rescues otherwise-unmatched
    reads with a 1 bp insertion or deletion inside the guide (see the module docstring);
    ignored with exact_match. Reads shorter than `min_len` or longer than `max_len` (None: no
    limit, the default) are skipped and reported as unmatched.

    A sample that can't be counted (no flanks learned, unreadable or malformed file) raises its
    error, unless `skip_failed`: then the other samples still run and that sample's entry in
    the returned list is the exception instead of a (counts, stats) tuple.

    `on_done(i, result)`, if given, is called in the calling process as each sample finishes,
    with its index in `paths` and its entry for the returned list, so progress can be shown
    before every sample is done.
    """
    paths = [os.fspath(p) for p in paths]
    if not paths:
        raise ValueError("No FASTQ files given.")
    seqs = [s.upper() for s in seqs]
    if not seqs:
        raise ValueError("No guide sequences given.")
    glen = len(seqs[0])
    if any(len(s) != glen for s in seqs):
        raise ValueError("Guide sequences must all be the same length.")
    if len(set(seqs)) != len(seqs):
        raise ValueError("Guide sequences must be unique.")
    for name, value, lo in (("flank_len", flank_len, 1), ("sample_size", sample_size, 1),
                            ("flank_hits", flank_hits, 1), ("min_len", min_len, 0),
                            ("max_len", max_len, max(1, min_len)), ("jobs", jobs, 1),
                            ("threads", threads, 1)):
        if value is not None and value < lo:
            raise ValueError(f"{name} must be at least {lo}, got {value}.")
    if max_len is None:
        max_len = _NO_MAX_LEN
    jobs = jobs or max(1, (os.cpu_count() or 2) // 2)
    threads = threads or max(1, (os.cpu_count() or 2) // max(1, min(jobs, len(paths))))
    work = [(p, seqs, glen, flank_len, sample_size, flank_hits, min_len, max_len, exact_match,
             threads, allow_indel) for p in paths]

    results = [None] * len(paths)

    def finish(i, get):
        results[i] = outcome(get, paths[i])
        if on_done is not None:
            on_done(i, results[i])

    def outcome(get, path):
        try:
            try:
                return get()
            except BrokenProcessPool as e:
                # the OS killed a worker (usually out of memory); every unfinished sample fails
                raise RuntimeError(
                    f"{path}: not counted because a worker process was killed, most likely "
                    "because it ran out of memory (each sample counted in parallel builds its own "
                    "guide lookup, roughly 0.75 GB for a genome-wide library). Run fewer samples "
                    "at once with a lower -j (--jobs).") from e
        except (RuntimeError, OSError) as e:
            if not skip_failed:
                raise
            return e

    if jobs == 1 or len(paths) == 1:
        for i, w in enumerate(work):
            finish(i, lambda w=w: count_sample(w))
        return results
    ex = ProcessPoolExecutor(max_workers=min(jobs, len(paths)), initializer=_ignore_sigint)
    try:
        index = {ex.submit(count_sample, w): i for i, w in enumerate(work)}
        pending = set(index)
        while pending:  # short waits, so Ctrl-C is noticed on every platform
            done, pending = wait(pending, timeout=0.5)
            for f in sorted(done, key=index.get):
                finish(index[f], f.result)
        return results
    except KeyboardInterrupt:
        # stop now: drop queued samples and end the running ones (the workers ignore Ctrl-C)
        procs = list((getattr(ex, "_processes", None) or {}).values())
        ex.shutdown(wait=False, cancel_futures=True)
        for proc in procs:
            proc.terminate()
        raise
    finally:
        ex.shutdown(wait=True)


def _ignore_sigint():
    """Pool worker initializer: Ctrl-C is handled by the main process, which stops the workers,
    instead of every worker printing its own traceback."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
