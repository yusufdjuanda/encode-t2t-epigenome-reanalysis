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

## Decisions that must be documented

1. Exact T2T-CHM13 assembly and annotation release.
2. ENCODE file-selection rules and matched controls.
3. Whether PCR duplicates are marked or removed.
4. Alignment flags and MAPQ policy.
5. Construction and application of the minimum-unique-k-mer masks.
6. Narrow- versus broad-peak settings for each target.
7. Coordinate conversion and comparison rules between references.
8. Deviations from the paper and reasons for each deviation.

## Completion criteria for the pilot

- The workflow runs from a clean clone using configuration only.
- All inputs have stable accessions and checksums.
- Every output records software versions and command lines.
- Failed jobs can be resumed safely.
- The final report contains alignment QC, peak QC, and reference-comparison
  summaries.
