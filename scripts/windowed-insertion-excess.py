import pysam
import numpy as np
from argparse import ArgumentParser


def main():
    args = get_args()
    samfile = pysam.AlignmentFile(args.bam)
    for line in open(args.bed):
        chrom, start, end = line.strip().split()
        insertion_excess = np.mean([get_insertion_excess(read.cigartuples)
                                    for read in samfile.fetch(chrom, int(start), int(end))])
        if not np.isnan(insertion_excess):
            print("{}:{}-{}\t{}".format(chrom, start, end, round(insertion_excess, ndigits=0)))


def get_insertion_excess(cigar):
    """Return excess insertions over deletions.

    Using pysam cigartuples
    sum all insertions and softclips (operation 1 and 4)
    minus sum of all deletions (operation 2)
    """
    return sum([l for o, l in cigar if o in [1, 4]]) - sum([l for o, l in cigar if o == 2])


def get_args():
    parser = ArgumentParser(description="count insertion excess per bin")
    parser.add_argument("bam", help="bam as created by jvarkit bamslicebed")
    parser.add_argument("bed", help="bed with windows")
    return parser.parse_args()


if __name__ == "__main__":
    main()
