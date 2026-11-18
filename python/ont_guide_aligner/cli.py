"""Command-line interface: ``ont-guide-aligner -l library.txt -o results/Screen *.fastq.gz``."""

import argparse
import os
import sys

from . import __version__
from ._core import close_guide_count
from .counter import LibraryError, count_fastqs, qc_stats, read_library, sample_name

DESCRIPTION = """\
Count CRISPR sgRNAs in Oxford Nanopore amplicon reads, one guide per read.

The vector flanks on each side of the guide are learned per sample, then found anywhere in
either orientation of each read, so guides are located despite indels, truncation or extra
sequence. A read goes to the single best-matching guide (up to 1 mismatch by default); reads
tied between two guides are dropped as ambiguous.

Library: tab- or comma-separated; guide_id, sequence, gene (all three required). With a
header row, columns are found by name in any order; without one, they are read in that order.

Outputs (prefix from --output):
  <output>.counts.txt            guide, gene, one count column per sample (MAGeCK-ready)
  <output>.alignment.stats.txt   per-sample totals, assignment rate, ambiguous/unmatched
                                 reads (unmatched includes reads skipped by
                                 --min/--max-read-len), learned flanks
  <output>.QC.stats.txt          per-sample library representation: zero-count guides,
                                 average/max count, Gini, P90/P10 skew ratio, % assigned
"""


# Below this % assigned, the run still writes its outputs but warns: that low a rate usually
# means a wrong library, the wrong FASTQ, or a --min/--max-read-len limit that skips most reads.
LOW_ASSIGNED_PCT = 50
# Warn about close guides only when more than this % of the library has another guide 1 base
# away (typical of base-editor and tiling libraries, not of knockout libraries).
CLOSE_GUIDE_PCT = 20
# Flanks learned from fewer exact guide matches than this, or that fewer than this share of the
# matches agree on, may be wrong; the run still finishes but warns.
MIN_FLANK_HITS = 20
MIN_FLANK_SHARE = 0.5


def int_at_least(lo):
    """argparse type: a whole number >= lo."""
    def parse(s):
        try:
            v = int(s)
        except ValueError:
            v = None
        if v is None or v < lo:
            raise argparse.ArgumentTypeError(f"must be a whole number of at least {lo}, got {s!r}")
        return v
    return parse


def os_reason(e):
    return e.strerror or str(e)


def warn_on_flanks(names, stats, flank_hits):
    """Warn when a sample's flanks were learned from few exact guide matches or without a clear
    majority, and when the samples learned different flanks from each other."""
    min_hits = min(MIN_FLANK_HITS, flank_hits)  # -H below the threshold stops learning sooner
    for name, st in zip(names, stats):
        if st["_flank_hits"] < min_hits:
            print(f"Warning: the vector flanks for {name} were learned from only "
                  f"{st['_flank_hits']} exact guide matches, so they may be wrong; the sample may "
                  "have too few reads to learn them reliably. Check its counts before using them.",
                  file=sys.stderr)
        if min(st["_up_frac"], st["_down_frac"]) < MIN_FLANK_SHARE:
            print(f"Warning: the learned vector flanks for {name} are shared by less than "
                  f"{100 * MIN_FLANK_SHARE:.0f}% of its exact guide matches (upstream "
                  f"{st['upstream_flank']}, downstream {st['downstream_flank']}). The reads may "
                  "come from more than one vector backbone (reads with the other flanks may go "
                  "uncounted) or be of poor quality.", file=sys.stderr)
    by_flanks = {}
    for name, st in zip(names, stats):
        by_flanks.setdefault(st["_flanks"], []).append(name)
    if len(by_flanks) > 1:
        groups = "; ".join(f"{up}|{down} in {', '.join(ns)}" for (up, down), ns in by_flanks.items())
        print(f"Warning: the samples learned different vector flanks: {groups}. Samples of one "
              "screen normally share one vector; a different flank may mean a different vector "
              "backbone or a mixed-up sample.", file=sys.stderr)


def build_parser():
    p = argparse.ArgumentParser(prog="ont-guide-aligner", description=DESCRIPTION,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastqs", nargs="+", help="Input FASTQ(.gz) files, one per sample.")
    p.add_argument("-l", "--library", required=True, help="Guide library file (id, sequence, gene).")
    p.add_argument("-o", "--output", required=True, help="Output prefix, e.g. results/MyScreen.")
    p.add_argument("-f", "--flank-len", type=int_at_least(1), default=6,
                   help="Length of the vector flank learned on each side of the guide (default 6).")
    p.add_argument("-N", "--sample-size", type=int_at_least(1), default=100000,
                   help="Maximum reads scanned to learn the vector flanks (default 100000).")
    p.add_argument("-H", "--flank-hits", type=int_at_least(1), default=2000,
                   help="Stop flank learning once this many exact guide matches have been seen "
                        "(default 2000).")
    p.add_argument("--exact-match", action="store_true",
                   help="Only count exact guide matches (default allows 1 mismatch). Use for "
                        "libraries whose guides differ by single bases, e.g. base-editor screens.")
    p.add_argument("--indel", action=argparse.BooleanOptionalAction, default=True,
                   help="Also count reads whose guide has a 1 bp insertion or deletion (common in "
                        "ONT reads, especially in homopolymers): when both flanks are found 1 bp "
                        "too far apart or too close, the bases between them must be one indel from "
                        "exactly one guide. Only reads that would otherwise be unmatched are "
                        "affected. On by default; --no-indel turns it off. Ignored with "
                        "--exact-match.")
    p.add_argument("-t", "--threads", type=int_at_least(1), default=None,
                   help="Threads per sample for a .gz made of many gzip members "
                        "(typical of ONT basecaller output). Default: CPUs divided by the number "
                        "of samples run at once.")
    p.add_argument("--min-read-len", type=int_at_least(0), default=0,
                   help="Skip reads shorter than this (default: no minimum). Skipped reads are "
                        "reported as unmatched.")
    p.add_argument("--max-read-len", type=int_at_least(1), default=None,
                   help="Skip reads longer than this (default: no maximum). Not normally needed: "
                        "a read of joined amplicons with two different guides goes to the guide "
                        "with more intact flanks (then exact over 1 mismatch), or is dropped as "
                        "ambiguous if they tie. Skipped reads are reported as unmatched.")
    p.add_argument("-j", "--jobs", type=int_at_least(1), default=max(1, (os.cpu_count() or 2) // 2),
                   help="Samples processed in parallel (default: half the CPUs).")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def write_outputs(prefix, ids, genes, names, results):
    """Write the three output files. Each is written to <name>.tmp first and all are moved into
    place only once all three are complete, so a failed write (disk full, ...) leaves earlier
    outputs untouched instead of a mix of old and new files."""
    targets = [f"{prefix}.counts.txt", f"{prefix}.alignment.stats.txt", f"{prefix}.QC.stats.txt"]
    tmps = [t + ".tmp" for t in targets]
    try:
        with open(tmps[0], "w", encoding="utf-8", newline="\n") as f:
            f.write("\t".join(["guide", "gene"] + names) + "\n")
            for i, (gid, gene) in enumerate(zip(ids, genes)):
                f.write("\t".join([gid, gene] + [str(r[0][i]) for r in results]) + "\n")
        write_table(tmps[1], names, [st for _, st in results])
        write_table(tmps[2], names, [qc_stats(c, st) for c, st in results])
        for tmp, target in zip(tmps, targets):
            os.replace(tmp, target)
    except BaseException:  # also Ctrl-C: don't leave .tmp files behind
        for tmp in tmps:
            try:
                os.remove(tmp)
            except OSError:
                pass
        raise


# Fixed decimal places in the stats files, trailing zeros kept (100.00, not 100.0).
DECIMALS = {
    "pct_assigned": 2,
    "Guide Representation (%)": 2,
    "Average Count": 2,
    "Gini Coefficient": 3,
    "Skew Ratio (P90/P10)": 2,
    "Percent Assigned (%)": 2,
}


def format_value(key, value):
    if key in DECIMALS and isinstance(value, (int, float)):
        return f"{value:.{DECIMALS[key]}f}"
    return str(value)


def write_table(path, names, rows):
    """One row per sample: a sample column, then each dict's keys as columns (except keys
    starting with "_", which carry values for the warnings only)."""
    keys = [k for k in rows[0] if not k.startswith("_")]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(["sample"] + keys) + "\n")
        for name, row in zip(names, rows):
            f.write("\t".join([name] + [format_value(k, row[k]) for k in keys]) + "\n")


def main(argv=None):
    p = build_parser()
    a = p.parse_args(argv)
    if a.max_read_len is not None and a.min_read_len > a.max_read_len:
        p.error(f"--min-read-len ({a.min_read_len}) is larger than --max-read-len ({a.max_read_len}).")
    if not os.path.basename(a.output):
        p.error(f"--output is a file name prefix such as results/MyScreen, not a folder "
                f"(got {a.output!r}).")
    missing = [f for f in a.fastqs if not os.path.isfile(f)]
    if missing:
        sys.exit(f"Error: input file not found: {missing[0]}"
                 + (f" (and {len(missing) - 1} more)" if len(missing) > 1 else ""))
    try:
        ids, seqs, genes, glen = read_library(a.library)
    except LibraryError as e:
        sys.exit(f"Error: {e}")
    except OSError as e:
        sys.exit(f"Error: can't read the library file {a.library}: {os_reason(e)}")
    except UnicodeDecodeError:
        sys.exit(f"Error: the library file {a.library} isn't a text file.")
    names = [sample_name(f) for f in a.fastqs]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        sys.exit(f"Error: two input files give the same sample name ({sorted(dup)[0]!r}).")

    out_dir = os.path.dirname(a.output)
    try:  # before counting, so a bad --output doesn't cost a long run
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
    except OSError as e:
        sys.exit(f"Error: can't create the output folder {out_dir}: {os_reason(e)}")

    print(f"Library: {len(ids)} guides of {glen} bp")
    print(f"Matching: {'exact only' if a.exact_match else 'up to 1 mismatch'}"
          + ("" if a.exact_match else f", 1 bp indel rescue {'on' if a.indel else 'off'}"))
    if not a.exact_match:
        close = close_guide_count(seqs)
        if close * 100 > CLOSE_GUIDE_PCT * len(seqs):
            sys.stdout.flush()
            print(f"Warning: {close:,} of {len(seqs):,} guides ({100 * close / len(seqs):.1f}%) "
                  "are only 1 base from another guide in the library, as in base-editor or "
                  "tiling libraries. With 1 mismatch allowed, a read with a sequencing error at "
                  "that base can be counted for the wrong guide. Use --exact-match if that "
                  "matters.", file=sys.stderr)

    def show(i, r):
        # one line per sample as it finishes (finishing order); outputs are written at the end
        if isinstance(r, Exception):
            line = f"{names[i]}: couldn't be counted (details at the end)"
        else:
            _, st = r
            line = (f"{names[i]}: {format_value('pct_assigned', st['pct_assigned'])}% assigned, "
                    f"flanks {st['upstream_flank']} | {st['downstream_flank']}")
        print(line, flush=True)

    try:
        results = count_fastqs(a.fastqs, seqs, exact_match=a.exact_match, flank_len=a.flank_len,
                               sample_size=a.sample_size, flank_hits=a.flank_hits,
                               min_len=a.min_read_len, max_len=a.max_read_len,
                               jobs=a.jobs, threads=a.threads, allow_indel=a.indel,
                               skip_failed=True, on_done=show)
    except (RuntimeError, OSError, ValueError) as e:
        sys.exit(f"Error: {e}")
    except KeyboardInterrupt:
        print("\nInterrupted; no output files were written.", file=sys.stderr)
        sys.exit(130)
    failed = [(n, r) for n, r in zip(names, results) if isinstance(r, Exception)]
    if len(failed) == len(names):
        sys.exit("Error: no sample could be counted.\n"
                 + "\n".join(f"  {n}: {e}" for n, e in failed))
    names, results = map(list, zip(*[(n, r) for n, r in zip(names, results)
                                     if not isinstance(r, Exception)]))

    try:
        write_outputs(a.output, ids, genes, names, results)
    except OSError as e:
        sys.exit(f"Error: can't write the output files {a.output}.*: {os_reason(e)}")
    out_dir, prefix = os.path.split(os.path.abspath(a.output))
    print(f"Wrote {prefix}.counts.txt, .alignment.stats.txt and .QC.stats.txt in {out_dir}")
    sys.stdout.flush()  # keep warnings after the summary when both streams go to one pipe
    limits = a.min_read_len > 0 or a.max_read_len is not None
    for name, (_, st) in zip(names, results):
        if st["pct_assigned"] < LOW_ASSIGNED_PCT:
            print(f"Warning: only {st['pct_assigned']:.2f}% of reads in {name} were assigned to a "
                  "guide. Check that the library matches this sample"
                  + (" and that --min-read-len / --max-read-len don't skip most reads "
                     "(skipped reads count as unmatched)" if limits else "")
                  + f". Learned flanks: {st['upstream_flank']} | {st['downstream_flank']}.",
                  file=sys.stderr)
    warn_on_flanks(names, [st for _, st in results], a.flank_hits)
    if failed:
        # outputs for the other samples are written; exit non-zero so pipelines notice
        sys.exit(f"Error: {len(failed)} of {len(failed) + len(names)} samples couldn't be counted "
                 "and are missing from the outputs:\n"
                 + "\n".join(f"  {n}: {e}" for n, e in failed))


if __name__ == "__main__":
    main()
