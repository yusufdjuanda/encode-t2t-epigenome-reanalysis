# ENCODE-T2T epigenome reanalysis

Reproducible reanalysis of ENCODE epigenomic sequencing data against the
T2T-CHM13 reference, inspired by Gershman et al., *Science* (2022),
"Epigenetic patterns in a complete human genome."

## Project status

The first milestone is a single ENCODE ChIP-seq experiment and its matched
input control (CTCF in C4-2B, `ENCSR460LGH`, with control `ENCSR585FQO`), run
end-to-end from raw FASTQ through alignment, filtering, peak calling and QC,
and compared with the paper's own results for the same sample.

Two analysis modes are planned:

- `paper_exact`: reproduce the published reference and software choices.
- `modern`: rerun the same samples with current reference/tool versions.

Results from the two modes must remain in separate output directories.
Everything below runs in `modern` mode; `paper_exact` is not set up yet.

### Done (pilot, `modern` mode)

- [x] Frozen accessions, FASTQ download and MD5 verification (12 FASTQs, 27 GB)
- [x] Alignment of all 6 library FASTQ pairs to T2T-CHM13 v2.0 (98.8-99.1%
      overall alignment rate for CTCF, 93.9-98.4% for the controls)
- [x] Merge into 4 per-replicate BAMs, with per-library manifests
- [x] The paper's methods and its authors' pipeline read and written up
      (`docs/paper-methods.md`), with the decisions and deviations recorded
      (`docs/replication-plan.md`)
- [x] Rules written for steps 03-06 (filter, deduplicate, unique-k-mer filter,
      pooled peaks, comparison), checked with a dry-run
- [x] The authors' KMC fork built (`build_kmc`): `kmc` and `kmc_genome_counts`
      in `<references>/tools/KMC-95c9190d4c/`

### To do

Nothing below has run yet.

1. **Build the unique-k-mer tracks** (`kmer_reference`, `count_kmers`,
   `make_kmer_track`): 7 k values, 23 sequences each, chrY excluded. Runtime,
   memory and output size are still unmeasured; the 256 GB requested for
   `wigToBigWig` is a guess. Try one k value before all seven.
2. **Run steps 03-05** for the pilot: filter and deduplicate each replicate,
   apply the k-mer filter, pool replicates per experiment, call peaks.
3. **Check the comparison** in
   `<results>/comparison/<reference>/paper_comparison.tsv`: read pairs kept at
   each step and the peak count, against the authors' C4-2B CTCF numbers
   (corrected for their double counting, see `docs/paper-methods.md`).
4. **Verify the k-mer filter** against the authors' `filter_by_unique_kmers.py`
   on a small BAM (deviation D10); the reimplementation is untested.
5. **Decide what `paper_exact` means here** now that CHM13 v2.0 is used for
   both modes (deviations D1, D2): at minimum the paper's bowtie2 flags, which
   need a second alignment run.
6. **Add the GRCh38 comparison** (phase 2): align the same libraries to
   GRCh38, run the same steps, lift peaks over to T2T and count the peaks
   found only in T2T. Needs a GRCh38 → CHM13 v2.0 chain file, and a check of
   which contigs the Karakoram GRCh38 bundle contains (deviation D6).
7. **Scale out** (phase 3) to more marks and biosamples from Table S1.
8. **Report**: alignment QC, peak QC and the reference comparison.

## Repository layout

```text
config/                     Version-controlled configuration
docs/                       Analysis decisions and replication notes
profiles/karakoram/         Snakemake profile for the Karakoram Slurm cluster
resources/accessions/       Frozen ENCODE experiment accessions (see below)
resources/                  Metadata only; large references are not committed
workflow/                   Snakemake workflow, rules, scripts, and environments
workflow/scripts/           One script per pipeline step, run by Snakemake rules
workflow/notebooks/         One exploration notebook per pipeline step
pixi.toml                   Pixi environments: `default` (Snakemake/pandas/Jupyter), `fetch` (ENCODEfetch)
```

Raw reads, genome indexes, BAM files, and results are intentionally excluded
from Git. Store them under `/storage/projects/$USER/` on Karakoram and use
`/scratch/$USER/$SLURM_JOB_ID` for temporary job data.

### Numbered pipeline steps

Each step of the pipeline gets a two-digit prefix shared by its script and
its notebook when a notebook is useful, e.g. `01_download_fastqs`:

```text
workflow/scripts/01_download_fastqs.py       Runs inside the Snakemake DAG
workflow/notebooks/01_download_fastqs.ipynb  Explores that step's output by hand
```

The prefix fixes the order libraries/data move through the pipeline and
makes it obvious which notebook to open to sanity-check a given step's
output. Numbers are assigned as steps are actually built, in the order laid
out in "Planned workflow" below — they are not reserved in advance.

Notebooks are plain `.ipynb` files, committed with outputs cleared (`jupyter
nbconvert --clear-output --inplace <file>.ipynb`) so re-running one doesn't
produce a noisy diff full of stale numbers. Open one in Jupyter (`pixi run
jupyter`, using the `default` environment's kernel — see Local setup) or VS
Code (point its Python interpreter at `.pixi/envs/default/bin/python`).
Every notebook resolves the repo root itself by walking up from its cwd to
find `pixi.toml` — Jupyter starts a notebook's kernel with the notebook's own
directory as cwd, not wherever the server was launched from, so paths can't
just assume cwd is the repo root. Per the reproducibility rules below,
promote anything you want to keep out of a notebook into a script or rule.

## Local setup

[pixi](https://pixi.sh) manages every tool this repo needs, in two isolated
environments defined by `pixi.toml`: `default` (Snakemake, pandas, Jupyter —
everything the workflow itself and its notebooks need) and `fetch`
(ENCODEfetch, kept separate so its dependencies never mix with Snakemake's).

```bash
pixi install
cp config/config.example.yaml config/config.yaml
pixi run snakemake --dry-run --printshellcmds
```

### Fetching ENCODE data

FASTQs are downloaded from the frozen ENCODE experiment accessions with
[ENCODEfetch](https://github.com/khan-lab/ENCODEfetch), which the
`download_fastqs` rule calls in pixi's `fetch` environment
(`pixi run -e fetch encodefetch ...`) so its dependencies stay isolated from
Snakemake's. The download rule also verifies MD5 checksums before downstream
alignment starts:

```bash
pixi run snakemake                                # download/verify FASTQs + align
pixi run snakemake data/fastq/download.verified   # download/verify FASTQs only
```

The accessions themselves are frozen in
`resources/accessions/pilot_accessions.txt` (one ENCODE experiment accession
per line; matched input controls are resolved automatically). FASTQ download
pulls real data (tens of GB even for the pilot) and belongs on Karakoram
under `/storage/projects`, not on a laptop.

### Aligning reads

Step 02 aligns every FASTQ pair in the samplesheet as its own job
(`align_library`: bowtie2 piped into `samtools sort`, then `samtools index`),
so libraries run in parallel on Slurm and a failure reruns only that library.
`download_fastqs` is a Snakemake checkpoint: the libraries to align are only
known once it has written the samplesheet. `alignment_manifest` then lists
every BAM with its sample, replicate and matched control for later steps.

```text
<results>/alignment/<reference>/<library>.sorted.bam(.bai)
<results>/alignment/<reference>/alignment_manifest.tsv
<logs>/alignment/<library>.bowtie2.log    bowtie2 alignment summary
<benchmarks>/alignment/<library>.tsv      runtime and memory of each job
```

Library IDs look like `ENCSR460LGH__case__rep1__ENCFF643KQW_ENCFF224VKL`
(experiment, case/control, replicate, read 1 and read 2 FASTQ accessions).
Request a single BAM to align just that library. The bowtie2 maximum fragment
length (`-X`) and the Lmod module versions are set under `alignment:` in the
config. `workflow/notebooks/02_align_fastqs.ipynb` walks through the step.

### Filtering, peaks and comparison with the paper

Steps 03-06 follow the authors' own pipeline (see `docs/paper-methods.md`,
with this project's deviations in `docs/replication-plan.md`):

| Step | Rules | Output |
|---|---|---|
| 03 Filter and deduplicate, per replicate | `filter_replicate` (samtools `-F 1804 -f 2 -q 2`), `deduplicate_replicate` (Picard, duplicates removed) | `<results>/filtering/<reference>/dedup/` |
| 04 Unique-k-mer filter, per replicate | `build_kmc`, `download_wigtobigwig`, `kmer_reference`, `count_kmers`, `make_kmer_track` (once per reference and k), `filter_unique_kmers` | tracks in `<references>/<reference>/kmer/`, BAMs in `<results>/filtering/<reference>/kmer/` |
| 05 Pooled peaks | `pool_experiment`, `call_peaks` (MACS3 `-f BAMPE`, case pooled vs control pooled) | `<results>/peaks/<reference>/macs3/` |
| 06 Comparison | `flagstat` per step, `compare_with_paper` | `<results>/comparison/<reference>/mapping_stats.tsv`, `paper_comparison.tsv` |

`paper_comparison.tsv` lists read pairs kept at each step (mapped, filtered,
dedup, kmer) and the peak count per experiment, next to the authors' numbers
from `resources/paper/T2T_Encode_Analysis`, corrected for their double
counting of multi-run libraries. The k-mer tracks are built once per reference
and reused by every mode.

## Karakoram setup

Clone the repository under project storage rather than the home directory:

```bash
cd /storage/projects/$USER
git clone <repository-url> encode-t2t-epigenome-reanalysis
cd encode-t2t-epigenome-reanalysis
cp config/config.example.yaml config/config.yaml
pixi install
```

Edit the copied config for Karakoram paths, including `paths.*`. Point
`reference.fasta` and `reference.bowtie2_index_prefix` either to the intended
shared bundle under `/storage/datasets/ref_genomes` or to project-storage
copies if the exact paper reference is not shared. Then test configuration
without starting work:

```bash
pixi run snakemake --profile profiles/karakoram --dry-run
```

Do not run alignment, indexing, or `download_fastqs` directly on the login
node. Use `--profile profiles/karakoram`, which submits work via Slurm.

## Planned workflow

1. Download raw FASTQs and verify checksums (ENCODEfetch).
2. Align paired-end reads to the configured T2T-CHM13 Bowtie2 index.
3. Retain primary alignments and mark/remove duplicates.
4. Apply the paper's unique-k-mer alignment filter.
5. Call peaks with the matched input using MACS2.
6. Generate alignment, library-complexity, and peak QC.
7. Annotate peaks unique to the T2T analysis and compare with GRCh38.

## Reproducibility rules

- Record accessions, checksums, reference versions, and command lines.
- Never overwrite `paper_exact` outputs with modernized results.
- Pin every workflow environment.
- Keep manual analysis out of the production workflow.
- Treat notebooks as exploration; promote final transformations into rules or
  scripts.

## Citation

Gershman A. et al. Epigenetic patterns in a complete human genome.
*Science*. 2022;376:eabj5089. doi:10.1126/science.abj5089.
