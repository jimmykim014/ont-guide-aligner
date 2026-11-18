# ont-guide-aligner

Count CRISPR sgRNAs in Oxford Nanopore (ONT) amplicon reads, **one guide per read**.

Tools built for short reads, such as guide-counter, look for the guide at fixed read offsets. On noisy ONT reads this can credit one read to several overlapping guides, or miss guides that indels have moved. ont-guide-aligner works in four steps:

1. It learns the vector sequence on each side of the guide from the first reads of each sample (6 bp each by default).
2. It searches for those flanks anywhere in both orientations of every read. This works for reads of any length, including truncated reads and reads with indels.
3. It looks up the window next to each flank in the library. Up to 1 mismatch is allowed by default; with `--exact-match`, only exact matches count. A 1-mismatch sequence that is 1 base from two different guides is never used; a read whose only hit is such a sequence is reported as ambiguous.
4. It gives each read to the single best guide. Two intact flanks beat one, and an exact match beats a 1-mismatch match. A read tied between two different guides is dropped as ambiguous; these are mostly two amplicons joined in one read.

Indel rescue (on by default; turn it off with `--no-indel`): a read that would otherwise be unmatched gets a second chance. This mostly recovers reads from guides with homopolymer runs, where ONT indels are most common.

## Install

```
pip install ont-guide-aligner
```

Requires Python 3.9 or later; there are no other dependencies. Prebuilt packages are available for Linux (x86_64 and ARM), macOS (Apple Silicon and Intel) and Windows (x64).

## Usage

```
ont-guide-aligner -l library.txt -o results/MyScreen sample1.fastq.gz sample2.fastq.gz ...
```

The library is tab- or comma-separated, one guide per row with a guide ID, a sequence and a gene. Guide IDs and sequences must be unique, and every guide must be the same length.

| Option                                  | Default         | Meaning                                                                                           |
| --------------------------------------- | --------------- | ------------------------------------------------------------------------------------------------- |
| `--exact-match`                       | off             | Count exact matches only (base-editor libraries)                                                  |
| `--indel` / `--no-indel`            | on              | Also count reads whose guide has a 1 bp insertion or deletion, common in homopolymers (see above); ignored with `--exact-match` |
| `-f/--flank-len`                      | 6               | Flank length learned on each side of the guide                                                    |
| `-H/--flank-hits`                     | 2000            | Stop flank learning after this many exact guide matches                                           |
| `-N/--sample-size`                    | 100000          | Scan at most this many reads to learn the flanks                                                  |
| `--min-read-len` / `--max-read-len` | none / none     | Skip reads outside this length range (counted as unmatched)                                      |
| `-j/--jobs`                           | half the CPUs   | Samples counted in parallel                                                                       |
| `-t/--threads`                        | CPUs ÷ samples run at once | Threads per sample for multi-member `.gz` files (typical ONT output)                            |

### Outputs

Each run writes three tab-separated files, named from the `-o` prefix.

**`<output>.counts.txt`**: one row per guide.

| Metric | Description |
| ------ | ----------- |
| `guide` | Guide ID from the library |
| `gene` | Gene from the library |
| one column per sample | Reads assigned to that guide in that sample |

**`<output>.alignment.stats.txt`**: how reads were assigned, one row per sample. `assigned_reads` + `ambiguous_dropped` + `unmatched` = `total_reads`.

| Metric | Description |
| ------ | ----------- |
| `sample` | Sample name (the FASTQ file name without its extensions) |
| `total_reads` | All reads in the file |
| `assigned_reads` | Reads assigned to a guide |
| `pct_assigned` | `assigned_reads` as a percentage of `total_reads` |
| `assigned_both_flanks` | Assigned reads with both flanks intact |
| `assigned_one_flank` | Assigned reads with only one flank intact (`assigned_both_flanks` + `assigned_one_flank` = `assigned_reads`) |
| `assigned_1_mismatch` | Assigned reads whose guide had 1 mismatch; already included in the two rows above |
| `assigned_indel` | Assigned reads recovered by indel rescue; already included in `assigned_both_flanks` (`NA` when the rescue is off) |
| `ambiguous_dropped` | Reads that fit more than one guide equally well |
| `unmatched` | Reads with no guide found, including reads skipped by `--min-read-len` / `--max-read-len` |
| `upstream_flank` / `downstream_flank` | The learned flanks, with how often each appeared next to exact guide matches while learning |

**`<output>.QC.stats.txt`**: library representation, one row per sample.

| Metric | Description |
| ------ | ----------- |
| `sample` | Sample name |
| `Total Guides` | Guides in the library |
| `Zero-Count Guides` | Guides with no reads |
| `Guide Representation (%)` | Guides with at least one read, as a percentage of `Total Guides` |
| `Total Counts` | Reads counted across all guides |
| `Average Count` | Mean reads per guide |
| `Max Count` | Reads for the most-counted guide |
| `Gini Coefficient` | Evenness of counts across guides: 0 is perfectly even, 1 is all reads on one guide (`NA` when no reads were counted) |
| `Skew Ratio (P90/P10)` | 90th / 10th percentile guide count (`NA` when the 10th percentile is 0) |
| `Total Reads` | All reads in the file |
| `Reads Assigned` | Reads assigned to a guide |
| `Percent Assigned (%)` | `Reads Assigned` as a percentage of `Total Reads` |

## Scope

This package is for ONT reads only. Each read must cover the guide plus some vector sequence on both sides. If no flanks can be learned, it stops with an error. That happens with short Illumina reads that end at the guide.

## Credits

The 1-mismatch lookup table and background decompression follow the approach of [guide-counter](https://github.com/fulcrumgenomics/guide-counter) (Fulcrum Genomics, MIT). No code was copied from it. The version whose approach informed this tool:

> Fennell T, Homer N. *fulcrumgenomics/guide-counter: v0.1.3*. Zenodo, 2022. [doi:10.5281/zenodo.6375792](https://doi.org/10.5281/zenodo.6375792)

The compiled core includes open-source Rust libraries; their licenses are in [THIRD_PARTY_LICENSES](https://github.com/jimmykim014/ont-guide-aligner/blob/main/THIRD_PARTY_LICENSES).

## Citation

If you use ont-guide-aligner, please cite it using the details in [CITATION.cff](https://github.com/jimmykim014/ont-guide-aligner/blob/main/CITATION.cff). GitHub also shows them under "Cite this repository".
