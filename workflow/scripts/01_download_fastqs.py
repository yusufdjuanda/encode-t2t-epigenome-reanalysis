"""Download ENCODE FASTQs from declared accessions and verify MD5 checksums."""

import csv
import hashlib
import subprocess
from pathlib import Path

accessions = Path(snakemake.input.accessions)
manifest_out = Path(snakemake.output.manifest)
metadata_out = Path(snakemake.output.metadata)
samplesheet_out = Path(snakemake.output.samplesheet)
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
with manifest_out.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("file_accession") and row.get("md5sum") and row.get("local_path"):
            expected[row["file_accession"]] = (row["md5sum"], Path(row["local_path"]))
        if row.get("file_accession_r2") and row.get("md5sum_r2") and row.get("local_path_r2"):
            expected[row["file_accession_r2"]] = (
                row["md5sum_r2"],
                Path(row["local_path_r2"]),
            )

if not expected:
    raise ValueError(f"No FASTQ paths/checksums found in manifest: {manifest_out}")

missing = [accession for accession, (_, path) in expected.items() if not path.exists()]
if missing:
    raise FileNotFoundError(f"Downloaded FASTQs missing for accessions: {sorted(missing)}")

mismatches = []
for accession, (want, path) in expected.items():
    got = md5sum(path)
    if got != want:
        mismatches.append((accession, want, got))

if mismatches:
    lines = "\n".join(f"{a}: expected {w}, got {g}" for a, w, g in mismatches)
    raise ValueError(f"Checksum mismatch for downloaded FASTQs:\n{lines}")

if not metadata_out.exists():
    raise FileNotFoundError(f"encodefetch did not write {metadata_out}")
if not samplesheet_out.exists():
    raise FileNotFoundError(f"encodefetch did not write {samplesheet_out}")

done_marker.write_text(f"Verified {len(expected)} FASTQ files against {manifest_out}.\n")
