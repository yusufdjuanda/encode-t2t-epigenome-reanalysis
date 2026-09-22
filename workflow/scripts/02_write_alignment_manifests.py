"""Write the per-library and per-replicate alignment manifests for downstream steps."""

import csv
from pathlib import Path

libraries = snakemake.params.libraries
replicates = snakemake.params.replicates
reference_name = snakemake.params.reference_name
library_bams = [Path(bam) for bam in snakemake.input.library_bams]
replicate_bams = [Path(bam) for bam in snakemake.input.replicate_bams]


def check_order(bams, ids, suffix):
    # The Snakefile lists BAMs in the same order as the IDs; check anyway so a
    # mismatch fails loudly instead of pairing a BAM with the wrong sample.
    if len(bams) != len(ids):
        raise ValueError(f"Got {len(bams)} BAMs for {len(ids)} entries")
    for bam, item in zip(bams, ids):
        if bam.name != f"{item}{suffix}":
            raise ValueError(f"BAM {bam} does not belong to {item}")


def write_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


check_order(library_bams, list(libraries), ".sorted.bam")
check_order(replicate_bams, list(replicates), ".sorted.bam")

write_tsv(
    snakemake.output.libraries,
    [
        "sample_id",
        "replicate_id",
        "sample",
        "group",
        "replicate",
        "control",
        "control_replicate",
        "antibody",
        "reference",
        "fastq_1",
        "fastq_2",
        "bam",
    ],
    [
        {
            "sample_id": library_id,
            "replicate_id": f"{row['sample']}__{row['group']}__rep{row['replicate']}",
            "sample": row["sample"],
            "group": row["group"],
            "replicate": row["replicate"],
            "control": row["control"],
            "control_replicate": row["control_replicate"],
            "antibody": row["antibody"],
            "reference": reference_name,
            "fastq_1": row["fastq_1"],
            "fastq_2": row["fastq_2"],
            "bam": bam,
        }
        for (library_id, row), bam in zip(libraries.items(), library_bams)
    ],
)

write_tsv(
    snakemake.output.replicates,
    [
        "replicate_id",
        "sample",
        "group",
        "replicate",
        "antibody",
        "control_replicate_id",
        "n_libraries",
        "library_ids",
        "reference",
        "bam",
    ],
    [
        {
            "replicate_id": replicate_id,
            "sample": entry["sample"],
            "group": entry["group"],
            "replicate": entry["replicate"],
            "antibody": entry["antibody"],
            "control_replicate_id": entry["control_replicate_id"],
            "n_libraries": len(entry["libraries"]),
            "library_ids": ",".join(entry["libraries"]),
            "reference": reference_name,
            "bam": bam,
        }
        for (replicate_id, entry), bam in zip(replicates.items(), replicate_bams)
    ],
)
