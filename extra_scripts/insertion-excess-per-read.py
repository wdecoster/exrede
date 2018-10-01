import pysam
from argparse import ArgumentParser


def main():
    args = get_args()
    samfile = pysam.AlignmentFile(args.bam)

    for read in samfile.fetch():
        insertion_excess = get_insertion_excess(read.cigartuples)
        print("{}\t{}".format(read.query_name, insertion_excess))


def get_insertion_excess(cigar):
    """Return excess insertions over deletions.

    Using pysam cigartuples
    sum all insertions (operation 1)
    minus sum of all deletions (operation 2)
    """
    return sum([l for o, l in cigar if o == 1]) - sum([l for o, l in cigar if o == 2])


def get_args():
    parser = ArgumentParser(description="report insertion excess per read")
    parser.add_argument("bam", help="bam")
    return parser.parse_args()


if __name__ == "__main__":
    main()
