import os
import sys

configfile: "config.yaml"


def get_samples(wildcards):
    return config["samples"][wildcards.sample]


def get_chromosomes(genome):
    '''
    Gets the chromosome identifiers from the fasta genome and bed annotation
    and returns the intersection of both
    '''
    fai = genome + ".fai"
    if not os.path.isfile(fai):
        sys.exit("Fasta index {} not found".format(fai))
    return list(set([i.split('\t')[0] for i in open(fai)]))


CHROMOSOMES = get_chromosomes(config["genome"], config["annotbed"])
