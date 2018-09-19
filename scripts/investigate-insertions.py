import pysam
import pandas as pd
import numpy as np
from argparse import ArgumentParser


class Locus(object):
    def __init__(self, locus):
        self.string = locus
        self.chrom = locus.split(':')[0]
        self.start, self.end = [int(i) for i in locus.split(':')[1].split('-')]
        self.length = self.end - self.start
        self.tup = (self.chrom, self.start, self.end)


class Insertion(object):
    def __init__(self, locus, excess):
        self.locus = Locus(locus)
        self.excess = excess


def main():
    args = get_args()
    df = pd.read_csv(args.bed, sep="\t", header=None, names=['locus', 'excess']) \
           .itertuples(index=False, name=None)
    bam = pysam.AlignmentFile(args.bam, "rb")
    next(df)  # skipping the header tuple
    for locusstring, excess in df:
        var = Insertion(locusstring, excess)
        median_coverage = get_coverage(args.bam, var.locus)
        reads_involved = get_reads(bam, var.locus)
        var.genotype = estimate_ploidy(len(reads_involved), median_coverage)
        get_inserted_sequence(reads_involved)


def get_args():
    parser = ArgumentParser(description="Investigate insertions in bam based on bed")
    parser.add_argument("bam", help="bam from which insertions were identified")
    parser.add_argument("bed", help="bed file in which intervals are identified with insertions")
    return parser.parse_args()


def get_coverage(bam, locus):
    """Return the median coverage in the interval"""
    return np.median([int(l.strip().split('\t')[2])
                      for l in pysam.samtools.depth("-r", locus.string, bam).split('\n') if l])


def get_reads(bam, locus):
    """Return all primary alignments contributing to the insertion event

    Parses cigartuples, entry[0] = operation, entry[1] is length
    Arbitrary cutoff of minimal insertion length = 50
    """
    reads_involved = []
    for read in bam.fetch(*locus.tup):
        if read.mapping_quality > 0:
            for entry in read.cigartuples:
                if entry[0] == 1 and entry[1] > 50:
                    reads_involved.append(read)
                    break
    return reads_involved


def estimate_ploidy(reads, coverage):
    ploidy = reads / coverage
    if 0.25 < ploidy < 0.75:
        return '0/1'
    elif ploidy >= 0.75:
        return '1/1'
    else:
        return './.'


def get_inserted_sequence(reads):


if __name__ == '__main__':
    main()
