"""Download the FASTQs resolved by 00_resolve_accessions.py via ENCODEfetch,
then verify every downloaded file's MD5 against the manifest ENCODEfetch
itself recorded. Intended to run on Karakoram via Slurm (large transfer),
not on a laptop or the login node.
"""

import csv
import hashlib
import subprocess
from pathlib import Path

accessions = Path(snakemake.input.accessions)
manifest_in = Path(snakemake.input.manifest)
done_marker = Path(snakemake.output.done)
outdir = done_marker.parent
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


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


expected = {}
with manifest_in.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("file_accession") and row.get("md5sum"):
            expected[row["file_accession"]] = row["md5sum"]
        if row.get("file_accession_r2") and row.get("md5sum_r2"):
            expected[row["file_accession_r2"]] = row["md5sum_r2"]

fastqs = {p.name.split(".")[0]: p for p in outdir.rglob("*.fastq.gz")}

missing = expected.keys() - fastqs.keys()
if missing:
    raise FileNotFoundError(f"Downloaded FASTQs missing for accessions: {sorted(missing)}")

mismatches = []
for accession, path in fastqs.items():
    want = expected.get(accession)
    if want is None:
        continue
    got = md5sum(path)
    if got != want:
        mismatches.append((accession, want, got))

if mismatches:
    lines = "\n".join(f"{a}: expected {w}, got {g}" for a, w, g in mismatches)
    raise ValueError(f"Checksum mismatch for downloaded FASTQs:\n{lines}")

done_marker.parent.mkdir(parents=True, exist_ok=True)
done_marker.write_text(f"Verified {len(fastqs)} FASTQ files against {manifest_in}.\n")
