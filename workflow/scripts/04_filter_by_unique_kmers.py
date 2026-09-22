"""Keep read pairs whose aligned reference span contains a genome-unique k-mer.

Re-implements bin/filter_by_unique_kmers.py from msauria/T2T_Encode_Analysis
(docs/paper-methods.md) with the same decision rule:

- Each k-mer track is a bigWig whose value at a position is the genome count of
  the k-mer starting there, so 1 means unique.
- A mate's span is reference_end - reference_start. Spans shorter than the
  smallest k fail. Otherwise k is the largest size strictly below the span,
  with the smallest size as the floor (sizes[searchsorted(sizes[1:-1], span)]),
  so a fully aligned 100 or 101 bp read uses k = 95.
- A mate passes if any k-mer starting inside its span, and ending inside it, is
  unique. Chromosomes absent from the tracks fail.
- A pair is kept if either mate passes.

Unlike the original, mates are paired by read name rather than by position in a
name-sorted file, so the input can stay coordinate-sorted. The output is in
pair order and needs sorting afterwards.
"""

import argparse
import bisect
import json

import pyBigWig
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bam", required=True, help="Paired-end BAM, primary proper pairs only")
    parser.add_argument("--out", required=True, help="Output BAM (unsorted)")
    parser.add_argument("--stats", required=True, help="Output JSON with pair counts")
    parser.add_argument("--sizes", required=True, help="Comma-separated k-mer sizes")
    parser.add_argument("--tracks", required=True, nargs="+", help="bigWig per size, same order")
    return parser.parse_args()


class UniqueKmerTracks:
    def __init__(self, sizes, paths):
        if len(sizes) != len(paths):
            raise ValueError(f"{len(sizes)} k-mer sizes but {len(paths)} tracks")
        order = sorted(range(len(sizes)), key=lambda i: sizes[i])
        self.sizes = [sizes[i] for i in order]
        self.tracks = [pyBigWig.open(paths[i]) for i in order]
        self.inner_sizes = self.sizes[1:-1]
        self.chromosomes = set(self.tracks[0].chroms())
        self.lookup_errors = 0

    def mate_passes(self, read):
        chrom = read.reference_name
        if chrom not in self.chromosomes:
            return False
        start = read.reference_start
        span = read.reference_end - start
        if span < self.sizes[0]:
            return False
        index = bisect.bisect_left(self.inner_sizes, span)
        k = self.sizes[index]
        # k-mer start positions whose k-mer lies entirely within the span.
        last_start = start + span - k + 1
        try:
            minimum = self.tracks[index].stats(chrom, start, last_start, type="min", exact=True)[0]
        except RuntimeError:
            self.lookup_errors += 1
            return False
        # Positions without a value are k-mers containing N; any stored value is a
        # count >= 1, so the minimum is 1 exactly when a unique k-mer is present.
        return minimum == 1

    def close(self):
        for track in self.tracks:
            track.close()


def main():
    args = parse_args()
    sizes = [int(size) for size in args.sizes.split(",")]
    tracks = UniqueKmerTracks(sizes, args.tracks)

    pairs_in = pairs_kept = 0
    pending = {}
    with pysam.AlignmentFile(args.bam, "rb") as bam, pysam.AlignmentFile(
        args.out, "wb", template=bam
    ) as out:
        for read in bam:
            mate = pending.pop(read.query_name, None)
            if mate is None:
                pending[read.query_name] = read
                continue
            pairs_in += 1
            if tracks.mate_passes(mate) or tracks.mate_passes(read):
                out.write(mate)
                out.write(read)
                pairs_kept += 1
    tracks.close()

    stats = {
        "pairs_in": pairs_in,
        "pairs_kept": pairs_kept,
        "unpaired_reads_dropped": len(pending),
        "track_lookup_errors": tracks.lookup_errors,
    }
    with open(args.stats, "w") as handle:
        json.dump(stats, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
