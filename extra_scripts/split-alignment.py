import mappy as mp
import gzip
from Bio import SeqIO
from argparse import ArgumentParser


class AlignedChunk(object):
    def __init__(self, start, strand):
        self.start = start
        self.strand = strand


def align(aligner, seq):
    '''
    Test if reads can get aligned to the lambda genome,
    if not: write to stdout
    '''
    for hit in aligner.map(seq):
        if hit.mapq > 0:
            return AlignedChunk(hit.r_st, hit.strand)


def make_chunks(seq, size=200):
    return [str(seq)[i:i + size] for i in range(0, len(str(seq)), size)]


def main():
    args = get_args()
    aligner = mp.Aligner(args.reference, preset="map-ont")
    for record in SeqIO.parse(gzip.open(args.reads, 'rt'), "fastq"):
        unique_chunks = [align(aligner, chunk) for chunk in make_chunks(record.seq)]
        aligned_chunks = [x for x in unique_chunks if x is not None]
        if aligned_chunks:
            ch0 = aligned_chunks[0]
            for ch in aligned_chunks:
                if abs(ch.start - ch0.start) > 1e6 and ch.strand != ch0.strand:
                    print(record.format("fastq"))
                    break
                ch0 = ch


def get_args():
    parser = ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("reads")
    return parser.parse_args()


if __name__ == '__main__':
    main()
