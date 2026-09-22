# Paper methods: ENCODE ChIP-seq reanalysis

What Gershman et al. (2022) did with ENCODE ChIP-seq, from two sources:

- **Methods text**: supplementary materials `science.abj5089_sm.pdf`, section 4.1
  "ENCODE Dynamic k-mer assisted mapping" (pages 36-37).
- **Authors' pipeline**: [msauria/T2T_Encode_Analysis](https://github.com/msauria/T2T_Encode_Analysis),
  commit `3837e8a` (2021-10-25; analysis code from `11bb0b0`, "Final analysis,
  first submission", 2021-05-26), GPL-3.0. A Snakemake workflow described as
  "Scripts for recapitulating the Encode reanalysis using the T2T chm13
  reference assembly". It builds its k-mer tracks with a KMC fork,
  [msauria/KMC](https://github.com/msauria/KMC) branch `kmer_mapping`.

Where the two disagree, the code is what produced the published numbers, so
this project follows the code and records the difference.

Other code the paper cites: [timplab/T2T-Epigenetics](https://github.com/timplab/T2T-Epigenetics)
(downstream figures only; no mapping code) and Zenodo
[10.5281/zenodo.6046354](https://doi.org/10.5281/zenodo.6046354) ("Code
repositories used for T2T Epigenetics"; not inspected, Zenodo returned
HTTP 504 on 2026-09-16).

## Pipeline, step by step

Paths `Snakefile:N` refer to the authors' Snakefile at `3837e8a`.

| # | Step | Paper / authors' code |
|---|---|---|
| 1 | Inputs | ENCODE ChIP-seq with ≥100 bp paired-end reads and a matched input control; accessions in their `data/encode_samples.txt` (= Table S1). |
| 2 | References | `chm13v1` = T2T-CHM13 **v1.0** from the UCSC T2T hub (`t2t-chm13-v1.0.fa.gz`, `Snakefile:237`); `GRCh38p13` = UCSC `hg38.fa.gz`, full file including alt, random and chrUn contigs (`Snakefile:224`). |
| 3 | Combine runs | FASTQs of the same ENCODE **library** concatenated before mapping (`rule concat`, `Snakefile:386`). |
| 4 | Align | Bowtie2 2.4.1, `--no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder`, default `-X 500`, output unsorted in read order (`Snakefile:454`). |
| 5 | Filter | samtools 1.10 `view -b -F 1804 -f 2 -q 2`, per library (`Snakefile:597`). |
| 6 | Duplicates | Picard MarkDuplicates 2.22.1 `VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=queryname REMOVE_DUPLICATES=true`, per library (`Snakefile:619`). |
| 7 | Pool | All libraries of an **experiment** (both replicates) merged with `samtools merge -n`, then coordinate-sorted (`rule concat_bam`, `Snakefile:641`). Pooling happens **before** the k-mer filter. |
| 8 | k-mer tracks | Per genome and k: `kmc -k{k} -m1024 -fm -ci2` (KMC fork, env pins `kmc 3.1.2rc1`), then `kmc_genome_counts` writes a wiggle of the genome count of the k-mer starting at each position, converted to bigWig (`Snakefile:476-553`). k = **50, 75, 80, 85, 90, 95, 100**. Chromosomes chr1-22 and chrX only. |
| 9 | k-mer filter | `bin/filter_by_unique_kmers.py` on the name-sorted pooled BAM (details below). |
| 10 | Tracks | deepTools 3.4.3 `bamCoverage -bs 1`; `bigwigCompare -bs 50` against the control. |
| 11 | Peaks | MACS2 2.2.7.1 `-f BAMPE -g 3.03e9` (chm13v1) or `2.79e9` (GRCh38p13); pooled experiment vs pooled control experiment. **Broad** peaks for H3K9me3, H3K27me3, H3K36me3; narrow for all other targets, including CTCF (`Snakefile:814-860`). |
| 12 | Liftover | UCSC `liftOver -minMatch=0.2 -bedPlus=3` with `hg38.t2t-chm13-v1.0.over.chain.gz`, peaks clipped to chromosome ends (`Snakefile:904`). |
| 13 | Comparison | `bedtools intersect -u -a lifted_GRCh38 -b chm13` counts overlapping peaks (`Snakefile:1194`, bedtools 2.26.0). CHM13-unique peaks: `bedtools subtract -A -a chm13 -b lifted_GRCh38` (timplab `encode/find_unique_peaks.sh`). |

### The k-mer filter in detail

Track values (`kmc_genome_counts.cpp`, KMC fork commit `95c9190`, the
"first submission" version): the value at a position is the genome count of
the k-mer starting there. `-ci2` leaves k-mers seen once out of the database,
and those positions are written as **1**. k-mers containing `N` are skipped.
KMC counts canonical k-mers by default, so both strands count.

Per read pair, in a name-sorted BAM (`filter_by_unique_kmers.py`):

1. For each mate, `span = reference_end - reference_start` (includes deletions,
   excludes insertions; bowtie2 end-to-end mode has no soft clipping).
2. `span < 50`: the mate fails.
3. Otherwise k = `sizes[searchsorted([75, 80, 85, 90, 95], span, side="left")]`,
   which is the **largest k strictly smaller than the span**, with 50 as the
   floor. A fully matched 100 or 101 bp read uses **k = 95**; the 100-mer track
   is never selected.
4. The mate passes if any k-mer start position in `[start, start + span - k]`
   has value 1 on that chromosome. A chromosome missing from the track (chrY,
   alt, random, chrUn, chrM) fails.
5. The **pair is kept if either mate passes**. Mates are taken as consecutive
   records, relying on name order.

## Where the methods text and the code disagree

| Methods text says | Code does |
|---|---|
| k-mer databases "for 50-100mers by multiples of five" | k = 50, 75, 80, 85, 90, 95, 100 |
| "If all 100bp of a read map then the 100mer database is used"; 99 bp → 95, 101 bp → 100 | Largest k **strictly** below the span: 100 or 101 bp → 95, 99 bp → 95 |
| "Mismatches only impact the k-mer size if they occur in the first or last positions" | Mismatches never change the span |
| KMC3 v3.1.1 | Custom KMC fork, conda env pins 3.1.2rc1 |
| "Alignments from replicates were then pooled" after the k-mer filter | Pooled after deduplication, before the k-mer filter (same result, since the filter is per pair) |
| MACS2 "with default parameters" | `-f BAMPE`, and `--broad` for H3K9me3, H3K27me3, H3K36me3 |
| Peaks "unique to CHM13 … using bedtools intersect" | `bedtools subtract -A` (equivalent) |
| "96 total sequencing libraries (table S9)" | Libraries are in Table S1; Table S9 lists CHM13 CDRs |

## Apparent double counting in the authors' pipeline

`EXPERIMENT_DICT[experiment]` receives one entry per **line** of
`encode_samples.txt`, and a library sequenced in two runs has two lines
(`Snakefile:18-27`). As a result:

- `rule concat_bam` passes the same library BAM twice to `samtools merge`, so
  the pooled `dedup_concat` and `dedup_kmer` BAMs contain every read pair of
  such libraries **twice**, after deduplication had already run.
- `rule get_mapping_stats` sums the per-library `mapped` and `filtered` counts
  once per line, counting those libraries twice.

This affects experiments with multi-run libraries, including the pilot
`ENCSR460LGH` (C4-2B CTCF: both libraries have two runs). Evidence: its
reported `mapped` count, 132,956,930 pairs, is larger than its total input,
98,787,303 pairs. The control `ENCSR585FQO` has one run per library and is not
affected. MACS2 in BAMPE mode keeps one copy of identical fragments by default
(`--keep-dup 1`), which likely limits the effect on peaks. Bigwig coverage and
Table S2 totals for affected experiments are inflated. This conclusion comes
from reading the code; it has not been re-run.

## Validation targets for the pilot

From the authors' `results/mapping_stats.txt` and `results/macs2_peak_overlap.txt`.
Counts are read pairs (flagstat reads ÷ 2). For CTCF, "as reported" is
double-counted as described above; "corrected" divides by 2.

**C4-2B CTCF, ENCSR460LGH** (total input 98,787,303 pairs)

| Step | chm13v1 as reported | chm13v1 corrected | GRCh38p13 as reported | GRCh38p13 corrected |
|---|---|---|---|---|
| mapped | 132,956,930 | 66,478,465 (67.3%) | 131,821,286 | 65,910,643 |
| filtered | 125,427,798 | 62,713,899 | 123,193,432 | 61,596,716 |
| dedup | 116,258,046 | 58,129,023 | 114,147,222 | 57,073,611 |
| kmer | 113,546,296 | 56,773,148 | 110,746,730 | 55,373,365 |

**C4-2B input control, ENCSR585FQO** (total input 64,646,039 pairs; not affected)

| Step | chm13v1 | GRCh38p13 |
|---|---|---|
| mapped | 44,519,560 (68.9%) | 43,918,208 |
| filtered | 40,900,944 | 40,097,661 |
| dedup | 40,167,515 | 39,351,078 |
| kmer | 39,388,485 | 38,266,308 |

**C4-2B CTCF peaks**: chm13v1 59,367; GRCh38p13 56,293; lifted to chm13v1
56,241; lifted peaks overlapping chm13v1 peaks 55,932.

Table S2 (per-target totals across all cell lines) is built from the same
statistics, so its CTCF and other multi-run rows include the double counting.
