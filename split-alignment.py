import mappy as mp
import gzip
from Bio import SeqIO
from argparse import ArgumentParser


def align(aligner, seq):
    '''
    Test if reads can get aligned to the lambda genome,
    if not: write to stdout
    '''
    for hit in aligner.map(seq):
        if hit.mapq > 0:
            return hit.r_st


def make_chunks(seq, size=200):
    return [str(seq)[i:i + size] for i in range(0, len(str(seq)), size)]


def main():
    args = get_args()
    aligner = mp.Aligner(args.reference, preset="map-ont")
    for record in SeqIO.parse(gzip.open(args.reads, 'rt'), "fastq"):
        unique_locations = [align(aligner, chunk) for chunk in make_chunks(record.seq)]
        aligned_locations = [x for x in unique_locations if x is not None]
        if aligned_locations:
            loc0 = aligned_locations[0]
            for loc in aligned_locations:
                if abs(loc - loc0) > 1e6:
                    print(record.format("fastq"))
                    break
                loc0 = loc


def get_args():
    parser = ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("reads")
    return parser.parse_args()


if __name__ == '__main__':
    main()
