# ENCODE-T2T epigenome reanalysis

Reproducible reanalysis of ENCODE epigenomic sequencing data against the
T2T-CHM13 reference, inspired by Gershman et al., *Science* (2022),
"Epigenetic patterns in a complete human genome."

## Project status

This repository is an initial scaffold. The first milestone is a single
ENCODE ChIP-seq experiment and its matched input control, run end-to-end from
raw FASTQ through alignment, filtering, peak calling, and QC.

Two analysis modes are planned:

- `paper_exact`: reproduce the published reference and software choices.
- `modern`: rerun the same samples with current reference/tool versions.

Results from the two modes must remain in separate output directories.

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
its notebook, e.g. `00_resolve_accessions`:

```text
workflow/scripts/00_resolve_accessions.py       Runs inside the Snakemake DAG
workflow/notebooks/00_resolve_accessions.ipynb  Explores that step's output by hand
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

ENCODE experiment accessions are resolved and downloaded with
[ENCODEfetch](https://github.com/khan-lab/ENCODEfetch), which the
`resolve_accessions`/`download_fastqs` rules call in pixi's `fetch`
environment (`pixi run -e fetch encodefetch ...`) so its dependencies stay
isolated from Snakemake's:

```bash
pixi run snakemake data/accessions/manifest.tsv   # resolve accessions -> metadata only
pixi run snakemake data/fastq/download.verified   # download FASTQs + verify checksums
```

The accessions themselves are frozen in
`resources/accessions/pilot_accessions.txt` (one ENCODE experiment accession
per line; matched input controls are resolved automatically). The download
step is not part of the default `snakemake` target — it pulls real FASTQ
data (tens of GB even for the pilot) and belongs on Karakoram under
`/storage/projects`, not on a laptop.

## Karakoram setup

Clone the repository under project storage rather than the home directory:

```bash
cd /storage/projects/$USER
git clone <repository-url> encode-t2t-epigenome-reanalysis
cd encode-t2t-epigenome-reanalysis
cp config/config.example.yaml config/config.yaml
pixi install
```

Edit the copied config for Karakoram paths, then test configuration without
starting work:

```bash
pixi run snakemake --profile profiles/karakoram --dry-run
```

Do not run alignment, indexing, or `download_fastqs` directly on the login
node — everything heavier than metadata resolution goes through
`--profile profiles/karakoram`, which submits via Slurm.

## Planned workflow

1. Resolve and freeze ENCODE experiment/file accessions (ENCODEfetch).
2. Download raw FASTQs and verify checksums (ENCODEfetch).
3. obtain the exact T2T-CHM13 reference used by the paper.
4. Build the Bowtie2 index.
5. Align paired-end reads.
6. retain primary alignments and mark/remove duplicates.
7. Apply the paper's unique-k-mer alignment filter.
8. Call peaks with the matched input using MACS2.
9. Generate alignment, library-complexity, and peak QC.
10. Annotate peaks unique to the T2T analysis and compare with GRCh38.

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
