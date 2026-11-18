import gzip
import os
import random
import re
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from ont_guide_aligner import cli
from ont_guide_aligner import counter as oga  # internal API, tested directly
from ont_guide_aligner.counter import revcomp

UP, DOWN = "ACACCG", "GTTTCA"


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def library(tmp_path):
    rng = random.Random(1)
    seqs = []
    while len(seqs) < 50:
        s = rand_seq(rng, 20)
        if s not in seqs:
            seqs.append(s)
    path = tmp_path / "lib.txt"
    path.write_text("guide_id\tsequence\tgene\n" +
                    "".join(f"g{i}\t{s}\tGENE{i // 5}\n" for i, s in enumerate(seqs)))
    return path, seqs


def write_fastq(path, reads, gz=False):
    text = "".join(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n" for i, r in enumerate(reads))
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        Path(path).write_text(text)


def amplicon(rng, guide, left=40, right=40):
    return rand_seq(rng, left - 6) + UP + guide + DOWN + rand_seq(rng, right - 6)


def mutate(s, pos):
    return s[:pos] + ("A" if s[pos] != "A" else "C") + s[pos + 1:]


def make_reads(seqs):
    """Reads with a known expected outcome, plus enough clean reads to learn the flanks."""
    rng = random.Random(2)
    expected = [0] * len(seqs)
    reads = []
    for k in range(400):  # clean reads, both orientations, varying guide position (ONT-like drift)
        j = k % len(seqs)
        r = amplicon(rng, seqs[j], left=30 + k % 7, right=30 + k % 5)
        reads.append(r if k % 2 else revcomp(r))
        expected[j] += 1
    reads.append(amplicon(rng, mutate(seqs[3], 10)))            # 1 mismatch -> g3
    expected[3] += 1
    reads.append(UP + seqs[4] + rand_seq(rng, 3))               # truncated right flank -> g4
    expected[4] += 1
    reads.append(amplicon(rng, seqs[5]) + amplicon(rng, seqs[6]))  # concatemer -> ambiguous
    reads.append(amplicon(rng, seqs[7]) + amplicon(rng, seqs[7]))  # same guide twice -> counts once
    expected[7] += 1
    reads.append(rand_seq(rng, 100))                            # unmatched
    reads.append("ACGT")                                         # too short -> unmatched
    return reads, expected


@pytest.mark.parametrize("gz", [False, True])
def test_counts_and_tallies(tmp_path, library, gz):
    _, seqs = library
    reads, expected = make_reads(seqs)
    fq = tmp_path / ("s.fastq.gz" if gz else "s.fastq")
    write_fastq(fq, reads, gz=gz)
    [(counts, st)] = oga.count_fastqs([fq], seqs, jobs=1, allow_indel=False)
    assert counts == expected
    assert st["total_reads"] == len(reads)
    assert "reads_outside_range" not in st
    assert st["ambiguous_dropped"] == 1
    assert st["unmatched"] == 2  # random read + the 4 bp read
    assert st["assigned_1_mismatch"] == 1
    assert st["assigned_one_flank"] == 1
    assert st["upstream_flank"].startswith(UP) and st["downstream_flank"].startswith(DOWN)
    # the column relationships documented in the README
    assert st["assigned_reads"] + st["ambiguous_dropped"] + st["unmatched"] == st["total_reads"]
    assert st["assigned_both_flanks"] + st["assigned_one_flank"] == st["assigned_reads"]


def test_read_length_limits(tmp_path, library):
    _, seqs = library
    reads, expected = make_reads(seqs)
    rng = random.Random(7)
    reads.append(rand_seq(rng, 300) + UP + seqs[9] + DOWN + rand_seq(rng, 300))  # 632 bp
    fq = tmp_path / "s.fastq"
    write_fastq(fq, reads)
    [(counts, st)] = oga.count_fastqs([fq], seqs, jobs=1)  # no limit by default: counted
    assert counts[9] == expected[9] + 1 and st["unmatched"] == 2
    [(counts, st)] = oga.count_fastqs([fq], seqs, jobs=1, max_len=500)  # skipped -> unmatched
    assert counts == expected and st["unmatched"] == 3


def test_flank_learning_respects_length_limits(tmp_path, library):
    _, seqs = library
    rng = random.Random(3)
    # long reads carry a different vector; only the short ones match the counted reads
    long_reads = [rand_seq(rng, 300) + "TTTTTT" + seqs[i % 50] + "AAAAAA" + rand_seq(rng, 300)
                  for i in range(1000)]
    fq = tmp_path / "s.fastq"
    write_fastq(fq, long_reads + make_reads(seqs)[0])
    [(_, st)] = oga.count_fastqs([fq], seqs, jobs=1, max_len=500)
    assert st["upstream_flank"].startswith(UP) and st["downstream_flank"].startswith(DOWN)


def test_cli_warns_on_low_assignment(tmp_path, library, capsys):
    lib, seqs = library
    rng = random.Random(4)
    reads = make_reads(seqs)[0]
    fq = tmp_path / "A.fastq"
    write_fastq(fq, reads + [rand_seq(rng, 150) for _ in range(len(reads))])  # ~50% junk
    cli.main(["-l", str(lib), "-o", str(tmp_path / "o"), "-j", "1", str(fq)])
    assert "Warning: only" in capsys.readouterr().err
    write_fastq(fq, reads)
    cli.main(["-l", str(lib), "-o", str(tmp_path / "o"), "-j", "1", str(fq)])
    assert "Warning" not in capsys.readouterr().err


def test_deletion_fitting_two_guides_is_ambiguous(tmp_path):
    rng = random.Random(11)
    core = rand_seq(rng, 20)
    seqs = [core[:10] + "A" + core[10:], core[:10] + "C" + core[10:]]  # 21 bp, 1 base apart
    while len(seqs) < 30:
        s = rand_seq(rng, 21)
        if s not in seqs:
            seqs.append(s)
    reads = [amplicon(rng, s) for s in seqs[2:] for _ in range(5)]  # exact reads to learn flanks
    reads.append(UP + core + DOWN)  # the 20 bp span is one deletion from both A and B
    fq = tmp_path / "s.fastq"
    write_fastq(fq, reads)
    [(counts, st)] = oga.count_fastqs([fq], seqs, jobs=1)
    assert counts[0] == counts[1] == 0
    assert st["ambiguous_dropped"] == 1 and st["unmatched"] == 0


def test_mismatch_fitting_two_guides_is_ambiguous(tmp_path):
    rng = random.Random(12)
    core = rand_seq(rng, 21)
    a = core[:5] + "A" + core[6:15] + "A" + core[16:]
    b = core[:5] + "C" + core[6:15] + "C" + core[16:]  # 2 bases from a
    seqs = [a, b]
    while len(seqs) < 30:
        s = rand_seq(rng, 21)
        if s not in seqs:
            seqs.append(s)
    reads = [amplicon(rng, s) for s in seqs[2:] for _ in range(5)]  # exact reads to learn flanks
    reads.append(UP + a[:15] + "C" + a[16:] + DOWN)  # 1 mismatch from both a and b
    fq = tmp_path / "s.fastq"
    write_fastq(fq, reads)
    for allow_indel in (True, False):
        [(counts, st)] = oga.count_fastqs([fq], seqs, jobs=1, allow_indel=allow_indel)
        assert counts[0] == counts[1] == 0 and sum(counts) == len(reads) - 1
        assert st["ambiguous_dropped"] == 1 and st["unmatched"] == 0


def test_gzip_detected_by_content(tmp_path, library):
    _, seqs = library
    reads, expected = make_reads(seqs)
    upper, misnamed = tmp_path / "S1.FASTQ.GZ", tmp_path / "S2.fastq"
    write_fastq(upper, reads, gz=True)
    write_fastq(misnamed, reads, gz=True)  # gzip content under a plain .fastq name
    for fq in (upper, misnamed):
        assert oga.count_fastqs([fq], seqs, jobs=1)[0][0] == expected
    assert oga.sample_name(str(upper)) == "S1"
    assert oga.sample_name("x/S2.fastq.bgz") == "S2"
    assert oga.sample_name("S3.Fq") == "S3"


def test_exact_match_drops_mismatch(tmp_path, library):
    _, seqs = library
    reads, expected = make_reads(seqs)
    fq = tmp_path / "s.fastq"
    write_fastq(fq, reads)
    [(counts, st)] = oga.count_fastqs([fq], seqs, exact_match=True, jobs=1)
    expected[3] -= 1
    assert counts == expected and st["assigned_1_mismatch"] == 0


def test_short_reads_without_flank_rejected(tmp_path, library):
    _, seqs = library
    fq = tmp_path / "illumina.fastq"
    write_fastq(fq, ["ACGTA" + seqs[i % 50] for i in range(200)])  # reads end at the guide
    with pytest.raises(RuntimeError, match="vector flanks"):
        oga.count_fastqs([fq], seqs, jobs=1)


def test_library_checks(tmp_path):
    p = tmp_path / "lib.txt"
    p.write_text("id\tseq\tgene\na\tACGTACGT\tX\na\tTTTTCCCC\tY\n")
    with pytest.raises(oga.LibraryError, match="IDs"):
        oga.read_library(p)
    p.write_text("id\tseq\tgene\na\tACGTACGT\tX\nb\tACGTACGT\tY\n")
    with pytest.raises(oga.LibraryError, match="sequences"):
        oga.read_library(p)
    p.write_text("id,seq,gene\na,ACGTACGT,X\nb,ACGTAC,Y\n")
    with pytest.raises(oga.LibraryError, match="same length"):
        oga.read_library(p)


def test_library_columns_and_header(tmp_path):
    p = tmp_path / "lib.txt"
    want = (["a", "b"], ["ACGTACGT", "TTTTCCCC"], ["X", "Y"], 8)
    p.write_text("a\tACGTACGT\tX\nb\tttttcccc\tY\n")  # no header: first guide kept
    assert oga.read_library(p) == want
    p.write_text("sequence,gene,guide_id\nACGTACGT,X,a\nTTTTCCCC,Y,b\n")  # columns by name
    assert oga.read_library(p) == want
    p.write_text("\ufeffGuide\tsgRNA Sequence\tGene\na\tACGTACGT\tX\nb\tTTTTCCCC\tY\n",
                 encoding="utf-8")
    assert oga.read_library(p) == want  # Excel BOM, other common names
    # the gene is required: no gene column, an empty gene cell, or a short headerless row
    p.write_text("id\tseq\na\tACGTACGT\nb\tTTTTCCCC\n")
    with pytest.raises(oga.LibraryError, match="no gene column"):
        oga.read_library(p)
    p.write_text("id\tseq\tgene\na\tACGTACGT\t\nb\tTTTTCCCC\tY\n")
    with pytest.raises(oga.LibraryError, match="line 2: missing the gene"):
        oga.read_library(p)
    p.write_text("a\tACGTACGT\tX\nb\tTTTTCCCC\n")
    with pytest.raises(oga.LibraryError, match="line 2: missing the gene"):
        oga.read_library(p)


def test_library_bad_rows(tmp_path):
    p = tmp_path / "lib.txt"
    p.write_text("id\tseq\tgene\na\tACGTACGT\tX\nguide_42\n")
    with pytest.raises(oga.LibraryError, match="line 3: missing the guide sequence"):
        oga.read_library(p)
    p.write_text("id\tseq\tgene\na\tACGTACGT\tX\nb\tACGNACGT\tY\n")
    with pytest.raises(oga.LibraryError, match="line 3: guide 'b'.*other than A, C, G, T \\(N\\)"):
        oga.read_library(p)
    p.write_text("sequence\tgene\nACGTACGT\tX\n")
    with pytest.raises(oga.LibraryError, match="no guide ID column"):
        oga.read_library(p)
    p.write_text("col1,col2,col3\nX,g1,ACGTACGT\n")  # unknown names, wrong order
    with pytest.raises(oga.LibraryError, match="taken in order"):
        oga.read_library(p)
    # no header and a bad first sequence: reported, not silently taken as a header
    p.write_text("g01,ACNT,GENE1\ng02,ACGT,GENE2\ng03,TTGG,GENE3\n")
    with pytest.raises(oga.LibraryError, match="line 1: guide 'g01'.*other than A, C, G, T"):
        oga.read_library(p)
    p.write_text("a,b,c\ng1,ACGT,X\ng2,TTGG,Y\n")  # unknown 1-letter names: still a header
    assert oga.read_library(p) == (["g1", "g2"], ["ACGT", "TTGG"], ["X", "Y"], 4)


def test_cli(tmp_path, library, capsys):
    lib, seqs = library
    reads, expected = make_reads(seqs)
    write_fastq(tmp_path / "A.fastq.gz", reads, gz=True)
    write_fastq(tmp_path / "B.fastq", reads)
    out = tmp_path / "res" / "Screen"
    cli.main(["-l", str(lib), "-o", str(out), "-j", "2",
              str(tmp_path / "A.fastq.gz"), str(tmp_path / "B.fastq")])
    lines = Path(f"{out}.counts.txt").read_text().splitlines()
    assert lines[0] == "guide\tgene\tA\tB"
    assert lines[1].split("\t") == ["g0", "GENE0", str(expected[0]), str(expected[0])]
    assert len(Path(f"{out}.alignment.stats.txt").read_text().splitlines()) == 3
    qc = [l.split("\t") for l in Path(f"{out}.QC.stats.txt").read_text().splitlines()]
    assert qc[0] == ["sample", "Total Guides", "Zero-Count Guides", "Guide Representation (%)",
                     "Total Counts", "Average Count", "Max Count", "Gini Coefficient",
                     "Skew Ratio (P90/P10)", "Total Reads", "Reads Assigned",
                     "Percent Assigned (%)"]
    row = dict(zip(qc[0], qc[1]))
    assert row["Total Guides"] == "50" and row["Total Reads"] == str(len(reads))
    assert row["Zero-Count Guides"] == str(expected.count(0))
    assert row["Total Counts"] == row["Reads Assigned"] == str(sum(expected))
    assert "Warning" not in capsys.readouterr().err  # normal samples: no flank warnings


def test_cli_warns_on_weak_flanks(tmp_path, library, capsys):
    lib, seqs = library
    rng = random.Random(8)
    reads, _ = make_reads(seqs)
    write_fastq(tmp_path / "A.fastq", reads)
    write_fastq(tmp_path / "Few.fastq", [amplicon(rng, seqs[i]) for i in range(10)])  # 10 matches
    out = tmp_path / "o"
    cli.main(["-l", str(lib), "-o", str(out), "-j", "1",
              str(tmp_path / "A.fastq"), str(tmp_path / "Few.fastq")])
    err = capsys.readouterr().err
    assert "flanks for Few were learned from only 10 exact guide matches" in err
    assert "for A " not in err and "different vector flanks" not in err
    # the values behind the warnings stay out of the stats files
    header = Path(f"{out}.alignment.stats.txt").read_text().splitlines()[0].split("\t")
    assert header == ["sample", "total_reads", "assigned_reads", "pct_assigned",
                      "assigned_both_flanks", "assigned_one_flank", "assigned_1_mismatch",
                      "assigned_indel", "ambiguous_dropped", "unmatched", "upstream_flank",
                      "downstream_flank"]
    # three vectors in equal parts: no flank is shared by half the matches
    vectors = [(UP, DOWN), ("TGTGGA", "CTGCAG"), ("GGATCC", "AAGCTT")]
    mixed = [rand_seq(rng, 30) + vectors[k % 3][0] + seqs[k % 50] + vectors[k % 3][1]
             + rand_seq(rng, 30) for k in range(300)]
    write_fastq(tmp_path / "Mixed.fastq", mixed)
    cli.main(["-l", str(lib), "-o", str(out), "-j", "1", str(tmp_path / "Mixed.fastq")])
    assert "flanks for Mixed are shared by less than 50%" in capsys.readouterr().err


def test_cli_warns_on_differing_flanks(tmp_path, library, capsys):
    lib, seqs = library
    rng = random.Random(9)
    write_fastq(tmp_path / "A.fastq", make_reads(seqs)[0])
    other = [rand_seq(rng, 30) + "TGTGGA" + seqs[k % 50] + "CTGCAG" + rand_seq(rng, 30)
             for k in range(200)]
    write_fastq(tmp_path / "B.fastq", [r if k % 2 else revcomp(r) for k, r in enumerate(other)])
    cli.main(["-l", str(lib), "-o", str(tmp_path / "o"), "-j", "1",
              str(tmp_path / "A.fastq"), str(tmp_path / "B.fastq")])
    err = capsys.readouterr().err
    assert ("the samples learned different vector flanks: ACACCG|GTTTCA in A; "
            "TGTGGA|CTGCAG in B.") in err
    assert "Warning: only" not in err and "learned from only" not in err


def test_failed_write_keeps_previous_outputs(tmp_path, library, monkeypatch):
    lib, seqs = library
    reads, _ = make_reads(seqs)
    write_fastq(tmp_path / "A.fastq", reads)
    out = tmp_path / "o"
    args = ["-l", str(lib), "-o", str(out), "-j", "1", str(tmp_path / "A.fastq")]
    cli.main(args)
    names = ["o.counts.txt", "o.alignment.stats.txt", "o.QC.stats.txt"]
    before = {n: (tmp_path / n).read_text() for n in names}
    write_fastq(tmp_path / "A.fastq", reads[:300])  # new counts, if they were written

    real_write_table = cli.write_table

    def failing_write_table(path, *rest):
        if "QC" in os.path.basename(path):  # the third file
            raise OSError(28, "No space left on device")
        return real_write_table(path, *rest)

    monkeypatch.setattr(cli, "write_table", failing_write_table)
    with pytest.raises(SystemExit) as e:
        cli.main(args)
    assert "can't write the output files" in str(e.value.code)
    assert "No space left on device" in str(e.value.code)
    assert {n: (tmp_path / n).read_text() for n in names} == before
    assert not list(tmp_path.glob("*.tmp"))


def test_library_other_header_names(tmp_path):
    p = tmp_path / "lib.txt"
    p.write_text("gene_id,UID,seq\nX,HGLibA_1,ACGTACGT\nY,HGLibA_2,TTTTCCCC\n",
                 encoding="utf-8")  # GeCKO: gene_id is the gene, UID the guide ID
    assert oga.read_library(p) == (["HGLibA_1", "HGLibA_2"], ["ACGTACGT", "TTTTCCCC"],
                                   ["X", "Y"], 8)
    p.write_text("Target Gene ID\tTarget Gene Symbol\tsgRNA ID\tsgRNA Target Sequence\n"
                 "1\tX\ta\tACGTACGT\n2\tY\tb\tTTTTCCCC\n", encoding="utf-8")  # Broad-style
    assert oga.read_library(p) == (["a", "b"], ["ACGTACGT", "TTTTCCCC"], ["X", "Y"], 8)


def test_library_utf16(tmp_path):
    """Excel's "Unicode Text" export: UTF-16 with a BOM, tab-separated."""
    p = tmp_path / "lib.txt"
    text = "guide_id\tsequence\tgene\r\na\tACGTACGT\tX\r\nb\tTTTTCCCC\tY\r\n"
    want = (["a", "b"], ["ACGTACGT", "TTTTCCCC"], ["X", "Y"], 8)
    for enc in ("utf-16-le", "utf-16-be"):
        bom = b"\xff\xfe" if enc == "utf-16-le" else b"\xfe\xff"
        p.write_bytes(bom + text.encode(enc))
        assert oga.read_library(p) == want


def test_qc_stats():
    counts = [0] + list(range(1, 10)) + [100]  # 11 guides
    st = {"total_reads": 200, "assigned_reads": 145, "pct_assigned": 72.5}
    qc = oga.qc_stats(counts, st)
    assert qc["Total Guides"] == 11 and qc["Zero-Count Guides"] == 1
    assert qc["Guide Representation (%)"] == 90.91
    assert qc["Total Counts"] == 145 and qc["Average Count"] == 13.18 and qc["Max Count"] == 100
    assert qc["Skew Ratio (P90/P10)"] == 9.0  # numpy.percentile: P90 = 9, P10 = 1
    assert qc["Total Reads"] == 200 and qc["Reads Assigned"] == 145
    assert oga.qc_stats([0, 0, 5, 9], st)["Skew Ratio (P90/P10)"] == "NA"
    assert oga.qc_stats([0, 0, 0], st)["Gini Coefficient"] == "NA"  # nothing counted


@pytest.mark.parametrize("args, message", [
    (["-f", "0"], "--flank-len: must be a whole number of at least 1, got '0'"),
    (["-j", "-1"], "--jobs: must be a whole number of at least 1"),
    (["-t", "0"], "--threads: must be a whole number of at least 1"),
    (["-N", "x"], "--sample-size: must be a whole number of at least 1, got 'x'"),
    (["--min-read-len", "300", "--max-read-len", "200"], "larger than --max-read-len"),
    (["-o", "results/"], "not a folder"),
])
def test_cli_rejects_bad_options(tmp_path, library, capsys, args, message):
    lib, _ = library
    with pytest.raises(SystemExit) as e:
        cli.main(["-l", str(lib), "-o", str(tmp_path / "out"), *args, str(tmp_path / "x.fastq")])
        # (a later -o overrides the first; the FASTQ needn't exist, options are checked first)
    assert e.value.code == 2 and message in capsys.readouterr().err


def test_cli_clean_errors(tmp_path, library):
    lib, seqs = library
    fq = tmp_path / "A.fastq"
    write_fastq(fq, make_reads(seqs)[0])
    blocker = tmp_path / "blocker"
    blocker.write_text("")

    def error(*args):
        with pytest.raises(SystemExit) as e:
            cli.main(list(args))
        return str(e.value.code)

    assert "can't read the library file" in error("-l", str(tmp_path / "nope.txt"), "-o",
                                                  str(tmp_path / "o"), str(fq))
    assert "input file not found" in error("-l", str(lib), "-o", str(tmp_path / "o"),
                                           str(tmp_path / "missing.fastq"))
    assert "can't create the output folder" in error("-l", str(lib), "-o",
                                                     str(blocker / "Screen"), str(fq))
    with pytest.raises(ValueError, match="flank_len must be at least 1"):
        oga.count_fastqs([fq], seqs, flank_len=0)


def test_cli_summary_line_decimals(tmp_path, library, capsys):
    lib, seqs = library
    write_fastq(tmp_path / "A.fastq", make_reads(seqs)[0])
    cli.main(["-l", str(lib), "-o", str(tmp_path / "o"), "-j", "1", str(tmp_path / "A.fastq")])
    line = next(l for l in capsys.readouterr().out.splitlines() if l.startswith("A: "))
    # always 2 decimals, as in the files; the other numbers are only in the stats files
    assert re.fullmatch(r"A: \d+\.\d\d% assigned, flanks \S+ \(\d+%\) \| \S+ \(\d+%\)", line), line


def test_cli_wrote_line(tmp_path, library, capsys, monkeypatch):
    lib, seqs = library
    write_fastq(tmp_path / "A.fastq", make_reads(seqs)[0])
    monkeypatch.chdir(tmp_path)
    cli.main(["-l", str(lib), "-o", "out/screen", "-j", "1", "A.fastq"])
    out = capsys.readouterr().out.splitlines()
    assert out[-1] == ("Wrote screen.counts.txt, .alignment.stats.txt and .QC.stats.txt in "
                       + str(tmp_path / "out"))
    assert (tmp_path / "out" / "screen.QC.stats.txt").exists()


def fastq_text(reads):
    return "".join(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n" for i, r in enumerate(reads))


def test_fastq_structure_checks(tmp_path, library):
    _, seqs = library
    reads, expected = make_reads(seqs)
    good = fastq_text(reads).splitlines(keepends=True)
    fq = tmp_path / "s.fastq"
    # blank lines between records and at the end are fine
    fq.write_text("".join(l + ("\n" if i % 4 == 3 else "") for i, l in enumerate(good)))
    assert oga.count_fastqs([fq], seqs, jobs=1)[0][0] == expected
    # a record missing its quality line, early (seen by flank learning) and late (Rust core)
    for drop in (11, len(good) - 9):
        fq.write_text("".join(good[:drop] + good[drop + 1:]))
        with pytest.raises((RuntimeError, OSError), match="not a valid FASTQ"):
            oga.count_fastqs([fq], seqs, jobs=1, flank_hits=5)
    # last quality line cut short: caught by the Rust core (flank learning stops early) and by
    # flank learning in Python when it scans the whole file
    for flank_hits in (5, 10**6):
        fq.write_text("".join(good[:-1]) + good[-1][:-5] + "\n")
        with pytest.raises((RuntimeError, OSError),
                           match="quality line is not the same length as its sequence"):
            oga.count_fastqs([fq], seqs, jobs=1, flank_hits=flank_hits)
    fq.write_text("".join(good[:-2]))  # ends in the middle of a read
    with pytest.raises(OSError, match="ends in the middle of a read"):
        oga.count_fastqs([fq], seqs, jobs=1)


@pytest.mark.parametrize("jobs", [1, 2])
def test_on_done_called_once_per_sample(tmp_path, library, jobs):
    _, seqs = library
    reads, expected = make_reads(seqs)
    paths = []
    for n in ("A", "B", "C"):
        write_fastq(tmp_path / f"{n}.fastq", reads)
        paths.append(tmp_path / f"{n}.fastq")
    rng = random.Random(5)
    write_fastq(tmp_path / "Empty.fastq", [rand_seq(rng, 100) for _ in range(50)])  # no guides
    paths.append(tmp_path / "Empty.fastq")
    seen = []
    res = oga.count_fastqs(paths, seqs, jobs=jobs, skip_failed=True,
                           on_done=lambda i, r: seen.append((i, r)))
    assert sorted(i for i, _ in seen) == [0, 1, 2, 3]
    assert all(res[i] is r for i, r in seen)  # the same objects that are returned
    assert all(res[i][0] == expected for i in range(3)) and isinstance(res[3], RuntimeError)


def test_failed_sample_is_skipped(tmp_path, library, capsys):
    lib, seqs = library
    reads, expected = make_reads(seqs)
    write_fastq(tmp_path / "A.fastq", reads)
    rng = random.Random(5)
    write_fastq(tmp_path / "Empty.fastq", [rand_seq(rng, 100) for _ in range(50)])  # no guides
    paths = [tmp_path / "A.fastq", tmp_path / "Empty.fastq"]
    with pytest.raises(RuntimeError, match="no exact guide matches"):
        oga.count_fastqs(paths, seqs, jobs=2)
    res = oga.count_fastqs(paths, seqs, jobs=2, skip_failed=True)
    assert res[0][0] == expected and isinstance(res[1], RuntimeError)
    out = tmp_path / "o"
    with pytest.raises(SystemExit) as e:
        cli.main(["-l", str(lib), "-o", str(out), "-j", "2", *map(str, paths)])
    assert "1 of 2 samples couldn't be counted" in str(e.value.code) and "Empty" in str(e.value.code)
    printed = capsys.readouterr().out.splitlines()
    assert "Empty: couldn't be counted (details at the end)" in printed
    assert any(l.startswith("A: ") and "% assigned" in l for l in printed)
    assert Path(f"{out}.counts.txt").read_text().splitlines()[0] == "guide\tgene\tA"
    assert len(Path(f"{out}.QC.stats.txt").read_text().splitlines()) == 2
    with pytest.raises(SystemExit) as e:  # every sample fails: no outputs, clear error
        cli.main(["-l", str(lib), "-o", str(tmp_path / "none"), str(paths[1])])
    assert "no sample could be counted" in str(e.value.code)
    assert not Path(f"{tmp_path / 'none'}.counts.txt").exists()


def test_parallel_gz_split_matches_single_thread(tmp_path, library, monkeypatch):
    """ONT .gz files are many gzip members; counting splits them across threads. That only
    happens for files of 64 MB or more, so lower the threshold to test it here."""
    from ont_guide_aligner import _core
    _, seqs = library
    reads, expected = make_reads(seqs)
    reads = reads * 30
    fq = tmp_path / "s.fastq.gz"
    text = fastq_text(reads).splitlines(keepends=True)
    with open(fq, "wb") as f:  # one gzip member per 200 lines (50 reads), like ONT output
        for i in range(0, len(text), 200):
            f.write(gzip.compress("".join(text[i:i + 200]).encode()))
    monkeypatch.setenv("ONT_GUIDE_ALIGNER_MIN_SPLIT_BYTES", "1")
    assert len(_core.gz_split_points(str(fq), 4)) == 4  # really split into 4 pieces
    [(c1, st1)] = oga.count_fastqs([fq], seqs, jobs=1, threads=1)
    [(c4, st4)] = oga.count_fastqs([fq], seqs, jobs=1, threads=4)
    assert c1 == c4 == [e * 30 for e in expected] and st1 == st4


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs a named pipe (POSIX)")
def test_ctrl_c_stops_counting_promptly(tmp_path, library):
    """Ctrl-C used to be ignored until the whole sample had been counted. The sample here is a
    named pipe fed for up to 10 s, so only an interrupt can end the count sooner. The Rust core
    is called directly (flanks given), so the pipe has a single reader."""
    from ont_guide_aligner import _core
    _, seqs = library
    chunk = fastq_text(make_reads(seqs)[0]).encode()
    fifo = tmp_path / "endless.fastq"
    os.mkfifo(fifo)
    (tmp_path / "chunk.fastq").write_bytes(chunk)
    # fed by a separate process, not a thread, which would need the GIL
    feeder = subprocess.Popen([sys.executable, "-c", (
        "import sys, time\n"
        "chunk = open(sys.argv[2], 'rb').read()\n"
        "end = time.monotonic() + 10\n"
        "try:\n"
        "    with open(sys.argv[1], 'wb') as f:\n"
        "        while time.monotonic() < end:\n"
        "            f.write(chunk)\n"
        "except BrokenPipeError:\n"
        "    pass\n"), str(fifo), str(tmp_path / "chunk.fastq")])
    timer = threading.Timer(0.5, signal.raise_signal, args=(signal.SIGINT,))
    t0 = time.monotonic()
    timer.start()
    try:
        with pytest.raises(KeyboardInterrupt):
            _core.count_file(str(fifo), seqs, len(seqs[0]), UP, DOWN, 0, 2**32 - 1, False, 1, True)
        assert time.monotonic() - t0 < 3  # without the fix: only after the 10 s feed ends
    finally:
        timer.cancel()
        feeder.kill()
        feeder.wait()


def _killed_worker(job):
    """Stands in for count_sample in a pool worker: dies the way an out-of-memory kill does."""
    os._exit(1)


def test_killed_worker_message(tmp_path, library, monkeypatch):
    _, seqs = library
    paths = [tmp_path / "A.fastq", tmp_path / "B.fastq"]
    for fq in paths:
        write_fastq(fq, make_reads(seqs)[0])
    monkeypatch.setattr(oga, "count_sample", _killed_worker)
    with pytest.raises(RuntimeError, match="worker process was killed.*out of memory.*lower -j"):
        oga.count_fastqs(paths, seqs, jobs=2)
    res = oga.count_fastqs(paths, seqs, jobs=2, skip_failed=True)
    assert all(isinstance(r, RuntimeError) and "out of memory" in str(r) for r in res)
    assert str(paths[0]) in str(res[0]) and str(paths[1]) in str(res[1])


def test_cli_warns_on_close_guides_only_above_threshold(tmp_path, library, capsys):
    lib, seqs = library  # 50 guides, none 1 base apart
    fq = tmp_path / "A.fastq"
    write_fastq(fq, make_reads(seqs)[0])

    def run(n_close, *extra):
        # n_close extra guides, each 1 base from one of g0..g(n_close-1)
        lib2 = tmp_path / f"lib_{n_close}.txt"
        lib2.write_text(Path(lib).read_text() + "".join(
            f"c{i}\t{mutate(seqs[i], 5)}\tCLOSE\n" for i in range(n_close)))
        cli.main(["-l", str(lib2), "-o", str(tmp_path / "o"), "-j", "1", *extra, str(fq)])
        return capsys.readouterr().err

    assert "1 base from another guide" not in run(0)
    assert "1 base from another guide" not in run(1)  # 2 of 51 guides (3.9%)
    assert "1 base from another guide" not in run(5)  # 10 of 55 (18.2%)
    assert "12 of 56 guides (21.4%) are only 1 base from another guide" in run(6)
    assert "1 base from another guide" not in run(6, "--exact-match")


def test_corrupt_gz_fails(tmp_path, library):
    _, seqs = library
    reads, _ = make_reads(seqs)
    fq = tmp_path / "s.fastq.gz"
    write_fastq(fq, reads, gz=True)
    data = fq.read_bytes()
    fq.write_bytes(data[: len(data) * 2 // 3])
    # the cut comes before 2000 guide matches, so flank learning (Python) meets it first
    with pytest.raises(RuntimeError, match=re.escape(str(fq)) + ": can't decompress"):
        oga.count_fastqs([fq], seqs, jobs=1)
    # flank learning stops after 5 matches: the Rust core meets it while counting
    with pytest.raises(OSError, match=re.escape(str(fq)) + ": failed to decompress"):
        oga.count_fastqs([fq], seqs, jobs=1, flank_hits=5)


def make_indel_reads(seqs):
    """make_reads plus reads whose guide carries a 1 bp deletion or insertion, with intact flanks."""
    reads, expected = make_reads(seqs)
    rng = random.Random(3)
    rescued = [0] * len(seqs)
    for k in range(60):
        j = (k * 7) % len(seqs)
        g = seqs[j]
        p = 1 + k % (len(g) - 2)
        g = g[:p] + g[p + 1:] if k % 2 else g[:p] + "T" + g[p:]
        r = rand_seq(rng, 34) + UP + g + DOWN + rand_seq(rng, 34)
        reads.append(r if k % 3 else revcomp(r))
        rescued[j] += 1
    return reads, expected, rescued


def test_indel_rescue(tmp_path, library):
    _, seqs = library
    reads, expected, rescued = make_indel_reads(seqs)
    fq = tmp_path / "s.fastq.gz"
    write_fastq(fq, reads, gz=True)
    [(off, st_off)] = oga.count_fastqs([fq], seqs, jobs=1, allow_indel=False)
    [(on, st_on)] = oga.count_fastqs([fq], seqs, jobs=1, allow_indel=True)
    # Without rescue, some indel reads can still match through a flank-adjacent window, so only
    # check that rescue adds exactly the indel reads that were otherwise unmatched.
    assert st_off["assigned_indel"] == "NA"
    assert st_on["assigned_indel"] == st_off["unmatched"] - st_on["unmatched"]
    assert all(b >= a for a, b in zip(off, on))  # rescue only adds reads
    assert on == [e + r for e, r in zip(expected, rescued)]
    assert st_on["ambiguous_dropped"] == st_off["ambiguous_dropped"]


def test_indel_ignored_with_exact_match(tmp_path, library):
    _, seqs = library
    reads, _, _ = make_indel_reads(seqs)
    fq = tmp_path / "s.fastq"
    write_fastq(fq, reads)
    a = oga.count_fastqs([fq], seqs, jobs=1, exact_match=True, allow_indel=False)
    b = oga.count_fastqs([fq], seqs, jobs=1, exact_match=True, allow_indel=True)
    assert a[0][0] == b[0][0] and b[0][1]["assigned_indel"] == "NA"


def test_cli_indel_default_on(tmp_path, library):
    lib, seqs = library
    reads, expected, rescued = make_indel_reads(seqs)
    write_fastq(tmp_path / "A.fastq", reads)
    for flag, want in (([], [e + r for e, r in zip(expected, rescued)]), (["--no-indel"], None)):
        out = tmp_path / ("on" if not flag else "off")
        cli.main(["-l", str(lib), "-o", str(out), "-j", "1", *flag, str(tmp_path / "A.fastq")])
        rows = [l.split("\t") for l in Path(f"{out}.counts.txt").read_text().splitlines()[1:]]
        stats = dict(zip(*[l.split("\t") for l in Path(f"{out}.alignment.stats.txt").read_text().splitlines()]))
        if want is not None:
            # some indel reads already match next to a flank without rescue, so assigned_indel <= 60
            assert [int(r[2]) for r in rows] == want and 0 < int(stats["assigned_indel"]) <= sum(rescued)
        else:
            assert stats["assigned_indel"] == "NA"
