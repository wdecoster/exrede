import pysam
from argparse import ArgumentParser


def main():
    args = get_args()
    samfile = pysam.AlignmentFile(args.bam)
    ins = [get_long_ins(read.cigartuples, read.query_sequence)
           for read in samfile.fetch() if read.query_sequence]
    head = 1
    for seq in ins:
        if seq:
            print(">ins{}".format(head))
            print(seq[0])
            head += 1


def get_insertion_excess(cigar):
    """Return excess insertions over deletions.

    Using pysam cigartuples
    sum all insertions (operation 1)
    minus sum of all deletions (operation 2)
    """
    return sum([l for o, l in cigar if o == 1]) - sum([l for o, l in cigar if o == 2])


def get_long_ins(cigartuples, sequence):
    cumsum = 0
    res = []
    for op, length in cigartuples:
        if op == 1 and length > 20:
            res.append(sequence[cumsum: cumsum + length])
        cumsum += length
    return res


def get_args():
    parser = ArgumentParser(description="")
    parser.add_argument("bam")
    return parser.parse_args()


if __name__ == "__main__":
    main()
