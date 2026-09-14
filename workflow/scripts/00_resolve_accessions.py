"""Resolve frozen ENCODE experiment accessions into file-level metadata via
ENCODEfetch (https://github.com/khan-lab/ENCODEfetch), without downloading
any FASTQs. Matched input-control experiments are resolved automatically.
"""

import subprocess
from pathlib import Path

accessions = Path(snakemake.input.accessions)
manifest_out = Path(snakemake.output.manifest)
metadata_out = Path(snakemake.output.metadata)
samplesheet_out = Path(snakemake.output.samplesheet)
outdir = manifest_out.parent
log_path = Path(snakemake.log[0])

assay_title = snakemake.params.assay_title
control_strategy = snakemake.params.control_strategy

outdir.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

with log_path.open("w") as log:
    subprocess.run(
        [
            "pixi",
            "run",
            "-e",
            "fetch",
            "encodefetch",
            "--accessions",
            str(accessions),
            "--assay-title",
            assay_title,
            "--file-type",
            "fastq",
            "--metadata-only",
            "--snakemake",
            "--control-strategy",
            control_strategy,
            "--no-progress",
            "--outdir",
            str(outdir),
        ],
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT,
    )

# ENCODEfetch names the samplesheet after --assay-title (e.g.
# snakemake_tf_chip-seq_samplesheet.csv); rename to the fixed path Snakemake
# expects.
(sheet,) = outdir.glob("snakemake_*_samplesheet.csv")
sheet.replace(samplesheet_out)

assert manifest_out.exists(), f"encodefetch did not write {manifest_out}"
assert metadata_out.exists(), f"encodefetch did not write {metadata_out}"
