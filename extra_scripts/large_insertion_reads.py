import pysam
from argparse import ArgumentParser


def main():
    args = get_args()
    samfile = pysam.AlignmentFile(args.bam)
    long_ins_file = pysam.AlignmentFile("long_ins.bam", "wb", template=samfile)
    for read in samfile.fetch():
        if get_insertion_excess(read.cigartuples) > 100:
            long_ins_file.write(read)


def get_insertion_excess(cigar):
    """Return excess insertions over deletions.

    Using pysam cigartuples
    sum all insertions (operation 1)
    minus sum of all deletions (operation 2)
    """
    return sum([l for o, l in cigar if o == 1]) - sum([l for o, l in cigar if o == 2])


def get_args():
    parser = ArgumentParser(description="count intertion excess per bin")
    parser.add_argument("bam")
    return parser.parse_args()


if __name__ == "__main__":
    main()
