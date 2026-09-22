"""Tabulate read pairs kept at each filtering step and peak counts, next to the paper's.

The authors' numbers come from resources/paper/T2T_Encode_Analysis (chm13v1).
Their mapping stats count a library once per sequencing run (deviation D4), so
they are divided by the runs per library where every library of the
experiment has the same number of runs.
"""

import csv
import json
from collections import defaultdict
from pathlib import Path

STEPS = ["mapped", "filtered", "dedup", "kmer"]

replicates = snakemake.params.replicates
experiments = snakemake.params.experiments
flagstat_keys = snakemake.params.flagstat_keys
peak_experiments = snakemake.params.peak_experiments
authors_dir = Path(snakemake.params.authors_dir)
reference_name = snakemake.params.reference_name


def read_pairs(path, total=False):
    counts = json.loads(Path(path).read_text())["QC-passed reads"]
    # Reads -> pairs. "primary" includes unaligned reads kept by bowtie2.
    return counts["primary" if total else "primary mapped"] // 2


def authors_mapping_stats():
    with (authors_dir / "mapping_stats.txt").open(newline="") as handle:
        rows = list(csv.reader(handle))
    genomes, header = rows[0], rows[1]
    table = {}
    for row in rows[2:]:
        if not row or not row[0]:
            continue
        entry = {"total": int(row[3])}
        for column, (genome, step) in enumerate(zip(genomes, header)):
            if genome == "chm13v1" and step in STEPS:
                entry[step] = int(row[column])
        table[row[0]] = entry
    return table


def authors_runs_per_library():
    """Runs per library for each experiment, or None when libraries differ."""
    runs = defaultdict(lambda: defaultdict(int))
    with (authors_dir / "encode_samples.txt").open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split(",")
            runs[fields[0]][fields[1]] += 1
    result = {}
    for experiment, libraries in runs.items():
        counts = set(libraries.values())
        result[experiment] = counts.pop() if len(counts) == 1 else None
    return result


def authors_labels():
    """Experiment accession -> the authors' track label, e.g. C4-2B_CTCF."""
    labels = {}
    with (authors_dir / "encode_samples.txt").open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split(",")
            cell_type = fields[6].replace(" ", "_").replace("(", "").replace(")", "")
            labels[fields[0]] = f"{cell_type}_{fields[5]}"
    return labels


def authors_peaks():
    with (authors_dir / "macs2_peak_overlap.txt").open(newline="") as handle:
        return {row["track"]: int(row["chm13v1"]) for row in csv.DictReader(handle, delimiter="\t")}


def percent(part, whole):
    return f"{100 * part / whole:.2f}" if part is not None and whole else ""


def main():
    flagstats = dict(zip((tuple(key) for key in flagstat_keys), snakemake.input.flagstats))

    # Per replicate, then summed per experiment.
    replicate_rows = []
    experiment_counts = defaultdict(lambda: defaultdict(int))
    for replicate_id, experiment_id in (
        (rep, exp) for exp, entry in experiments.items() for rep in entry["replicates"]
    ):
        counts = {"total": read_pairs(flagstats[("mapped", replicate_id)], total=True)}
        for step in STEPS:
            counts[step] = read_pairs(flagstats[(step, replicate_id)])
        replicate_rows.append({"replicate_id": replicate_id, "experiment_id": experiment_id, **counts})
        for key, value in counts.items():
            experiment_counts[experiment_id][key] += value

    fields = ["replicate_id", "experiment_id", "total", *STEPS]
    with open(snakemake.output.mapping_stats, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(replicate_rows)
        for experiment_id, counts in experiment_counts.items():
            writer.writerow({"replicate_id": "all", "experiment_id": experiment_id, **counts})

    peak_counts = {}
    for experiment_id, path in zip(peak_experiments, snakemake.input.peaks):
        with open(path) as handle:
            peak_counts[experiment_id] = sum(1 for line in handle if line.strip())

    stats = authors_mapping_stats()
    runs = authors_runs_per_library()
    labels = authors_labels()
    paper_peaks = authors_peaks()

    fields = [
        "experiment_id", "accession", "antibody", "measure",
        "this_project", "this_project_pct_of_total",
        "paper_reported", "paper_runs_per_library", "paper_corrected", "paper_pct_of_total",
        "this_project_vs_paper_pct",
    ]
    with open(snakemake.output.comparison, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for experiment_id, entry in experiments.items():
            accession = entry["sample"]
            ours = experiment_counts[experiment_id]
            theirs = stats.get(accession)
            divisor = runs.get(accession)
            for measure in ["total", *STEPS]:
                reported = theirs.get(measure) if theirs else None
                if reported is None:
                    corrected = None
                elif measure == "total":
                    # Their total comes from ENCODE metadata, not from BAMs.
                    corrected = reported
                else:
                    corrected = reported // divisor if divisor else None
                writer.writerow({
                    "experiment_id": experiment_id,
                    "accession": accession,
                    "antibody": entry["antibody"],
                    "measure": f"{measure}_pairs",
                    "this_project": ours[measure],
                    "this_project_pct_of_total": percent(ours[measure], ours["total"]),
                    "paper_reported": "" if reported is None else reported,
                    "paper_runs_per_library": "" if measure == "total" or divisor is None else divisor,
                    "paper_corrected": "" if corrected is None else corrected,
                    "paper_pct_of_total": percent(corrected, theirs["total"] if theirs else None),
                    "this_project_vs_paper_pct": (
                        f"{100 * (ours[measure] - corrected) / corrected:+.2f}" if corrected else ""
                    ),
                })
            if experiment_id in peak_counts:
                paper = paper_peaks.get(labels.get(accession, ""))
                writer.writerow({
                    "experiment_id": experiment_id,
                    "accession": accession,
                    "antibody": entry["antibody"],
                    "measure": f"peaks_{reference_name}_vs_chm13v1",
                    "this_project": peak_counts[experiment_id],
                    "paper_reported": "" if paper is None else paper,
                    "paper_corrected": "" if paper is None else paper,
                    "this_project_vs_paper_pct": (
                        f"{100 * (peak_counts[experiment_id] - paper) / paper:+.2f}" if paper else ""
                    ),
                })


main()
