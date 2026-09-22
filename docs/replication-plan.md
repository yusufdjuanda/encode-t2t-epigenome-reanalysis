# Replication plan

## Primary question

Can raw ENCODE ChIP-seq data, when realigned to T2T-CHM13 and filtered using
the published unique-k-mer strategy, reproduce the reported increase in
detectable epigenomic peaks and reveal signal in regions unresolved by GRCh38?

## Phase 1: pilot

- Select one ENCODE biosample.
- Select one histone mark and its matched input control.
- Freeze experiment and file accessions.
- Run FASTQ download, checksum verification, alignment, duplicate handling,
  unique-k-mer filtering, MACS2 peak calling, and QC.
- Inspect read retention and peak counts before scaling.

## Phase 2: reference comparison

- Process the same raw libraries against GRCh38 and T2T-CHM13.
- Use identical filtering and peak-calling parameters where valid.
- Compare mapped reads, uniquely retained reads, duplicate rates, peak counts,
  and genomic annotations.

## Phase 3: scale-out

- Add marks and biosamples from the publication.
- Use Slurm job arrays or Snakemake's Slurm executor.
- Stage high-I/O temporary work under `/scratch/$USER/$SLURM_JOB_ID`.
- Copy final BAMs, peaks, QC, and provenance back to `/storage/projects`.

## Decisions

What the paper did is recorded in [paper-methods.md](paper-methods.md), from
the supplementary methods and the authors' own pipeline
(msauria/T2T_Encode_Analysis). Where the two disagree, the authors' code is
treated as the reference.

| # | Decision | Status |
|---|---|---|
| 1 | T2T-CHM13 release | **Decided (2026-09-16): T2T-CHM13 v2.0 for all runs.** The paper used v1.0. See deviation D1. |
| 2 | ENCODE file selection and matched controls | Paper: ≥100 bp paired-end, matched input control, runs of one library combined (Table S1). This project: ENCODEfetch with `--control-strategy best`. |
| 3 | PCR duplicates | Paper: removed, Picard `REMOVE_DUPLICATES=true`, per library. |
| 4 | Alignment flags and MAPQ policy | Paper: bowtie2 `--no-discordant --no-mixed --very-sensitive --no-unal --omit-sec-seq --xeq --reorder` (default `-X 500`); samtools `-F 1804 -f 2 -q 2`. |
| 5 | Unique-k-mer filter | Paper: k = 50, 75-100 by 5; largest k below the aligned reference span; pair kept if either mate has a unique k-mer; tracks cover chr1-22 and X. This project: same, tracks built from v2.0 without chrY (D1). |
| 6 | Narrow vs broad peaks | Paper: MACS2 `-f BAMPE`, broad for H3K9me3, H3K27me3, H3K36me3, narrow otherwise (CTCF narrow). This project: same, with MACS3 (D7) and genome size 2,814,334,875 (D9). |
| 7 | Coordinate conversion and comparison | Paper: liftOver GRCh38 peaks with `-minMatch=0.2`; overlap with `bedtools intersect -u`; CHM13-unique with `bedtools subtract -A`. Chain for v2.0: open (D1). |
| 8 | Deviations from the paper | Listed below. |

## Deviations from the paper

| ID | Paper (authors' code) | This project | Reason | Consequences / to do |
|---|---|---|---|---|
| D1 | T2T-CHM13 **v1.0** (UCSC hub FASTA, `chr1`… names, no chrY) | T2T-CHM13 **v2.0**, Karakoram shared bundle `/storage/datasets/ref_genomes/t2t_chm13v2` (RefSeq names `NC_060925.1`…, includes the HG002 chrY) | Decided 2026-09-16: current release, already indexed on Karakoram. | Mapped counts on chrY and any sequence changed since v1.0 won't match the paper exactly. k-mer tracks must be built from v2.0. The v1.0 liftover chain can't be reused; a GRCh38 → v2.0 chain is needed. **Decided 2026-09-16: k-mer tracks exclude chrY**, matching the paper's chr1-22 and X (reads on chrY then fail the k-mer filter). Alignment still uses the full v2.0 index including chrY. Open: MACS2 genome size for v2.0 (paper used 3.03e9 for v1.0). |
| D2 | Bowtie2 2.4.1 with the flags in decision 4, `-X 500`, output in read order | `modern` mode: bowtie2 2.5.5 module, default sensitivity, `-X 2000`, discordant and unaligned reads kept, coordinate-sorted, read groups added | `-X 500` with `--no-discordant` drops pairs with fragments over 500 bp (~30% in a 20,000-pair test of ENCFF643KQW/ENCFF224VKL). | Current BAMs are not comparable at the "mapped" step. A `paper_exact` alignment with the paper's flags is still to build. |
| D3 | FASTQs of a library concatenated, aligned once; libraries pooled per experiment | Each run aligned separately, then merged per biological replicate | Parallel per-run jobs; per-run QC. | Same reads per replicate when one replicate = one ENCODE library (true for the pilot). Pooling per experiment for peak calling is still to build. |
| D4 | Multi-run libraries counted twice in pooled BAMs and mapping stats (see paper-methods.md) | Not reproduced: the runs of a library are merged once, so each read pair is counted once | **Decided 2026-09-16.** Appears to be a bug in the authors' pipeline; the paper's intent is to combine a library's runs. | Compare pilot CTCF results against the corrected (÷2) numbers in paper-methods.md. |
| D5 | Custom KMC fork (`kmc_genome_counts`), conda KMC 3.1.2rc1 | Not built yet | — | Decide whether to build the fork or reimplement the per-position unique-k-mer track. |
| D6 | GRCh38p13 = UCSC `hg38.fa.gz` including alt, random and chrUn contigs | Not built yet; Karakoram bundle `/storage/datasets/ref_genomes/grch38` not yet checked | — | Check the bundle's contig set; alt contigs lower MAPQ and change the `-q 2` filter. |
| D7 | MACS2 2.2.7.1 | MACS3 3.0.4 (Karakoram module), same arguments (`callpeak -f BAMPE`, `--broad` for broad marks) | MACS2 is not available as a module; `modern` mode. | Peak counts may differ slightly between versions. |
| D8 | Picard 2.22.1 MarkDuplicates on name-ordered BAMs (`ASSUME_SORT_ORDER=queryname`), per library | Picard 3.4.0 on coordinate-sorted BAMs, per replicate, `--REMOVE_DUPLICATES true` | Alignment output is coordinate-sorted; `modern` mode. | Coordinate mode can mark slightly different pairs. Read groups give every run of a replicate the same LB, so duplicates are found across runs as in the paper. |
| D9 | MACS2 `-g 3.03e9` (chm13v1) | `-g 2814334875` | Decided 2026-09-16: deepTools' effective genome size for T2T-CHM13 v2 at 100 bp reads, meant for MAPQ-filtered data. | Lower genome-wide background than the paper's; local backgrounds from the control usually dominate. |
| D10 | Unique-k-mer filter: `bin/filter_by_unique_kmers.py` on a name-sorted BAM, `pyBigWig.values()`; samtools 1.10 | `workflow/scripts/04_filter_by_unique_kmers.py`: same decision rule, mates paired by read name on the coordinate-sorted BAM, `pyBigWig.stats(type="min")`; per replicate, pooled afterwards; samtools 1.23.1 | Avoids a name sort; `modern` mode. | Should give identical pairs; not yet verified against the original script. |
| D11 | k-mer tracks from `kmc_genome_counts` in the authors' KMC fork, conda KMC 3.1.2rc1 build | Same fork at commit `95c9190d4c`, built from source with the system gcc 13 (rule `build_kmc`) | Only source of `kmc_genome_counts`. | The 2021 code has not been compiled with gcc 13 yet; the build may need fixes. |

## Completion criteria for the pilot

- The workflow runs from a clean clone using configuration only.
- All inputs have stable accessions and checksums.
- Every output records software versions and command lines.
- Failed jobs can be resumed safely.
- The final report contains alignment QC, peak QC, and reference-comparison
  summaries.
