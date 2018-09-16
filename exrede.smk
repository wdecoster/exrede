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


CHROMOSOMES = get_chromosomes(config["genome"])

rule all:
    input:
        expand("success-{sample}.txt", sample=config["samples"])


rule makewindows:
    input:
        fasta = config["genome"],
        fai = config["genome"] + ".fai"
    params:
        windowsize = config["windowsize"]
    output:
        "utils/windows.bed"
    log:
        "logs/makewindows.log"
    shell:
        """
        bedtools makewindows -g <(cut -f1-2 {input.fai}) -w {params.windowsize} | \
         awk -F $'\\t' '{{printf("%s\\t%s\\t%d\\n",$1,$2,int($3)-1);}}' > {output} 2> {log}
        """

rule samtools_split:
    input:
        bam = get_samples,
    output:
        temp("alignment/{sample}-{chromosome}.bam")
    params:
        chrom = "{chromosome}"
    log:
        "logs/samtools_split/{sample}-{chromosome}.log"
    shell:
        "samtools view {input.bam} {params.chrom} -o {output} 2> {log}"


rule samtools_index:
    input:
        "{dir}/{sample}.bam"
    output:
        "{dir}/{sample}.bam.bai"
    threads: 4
    log:
        "logs/samtools_index/{sample}.log"
    shell:
        "samtools index -@ {threads} {input} 2> {log}"


rule bamslicebed:
    input:
        bed = "utils/windows.bed",
        bam = "alignment/{sample}-{chromosome}.bam",
        bai = "alignment/{sample}-{chromosome}.bam.bai",
    output:
        temp("alignment_sliced/{sample}-{chromosome}.bam"),
    log:
        "logs/bamslicedbed/{sample}-{chromosome}.log"
    params:
        chrom = "{chromosome}"
    threads: 8
    shell:
        """
        bamslicebed -B <(grep ^{params.chrom} {input.bed}) {input.bam} | \
         samtools sort -@ {threads} -o {output} 2> {log}
        """

rule windowed_insertion_excess:
    input:
        bam = "alignment_sliced/{sample}-{chromosome}.bam",
        bed = "utils/windows.bed",
    output:
        "insertion_excess/{sample}-{chromosome}.ie"
    log:
        "logs/windowed-insertion-excess/{sample}-{chromosome}.log"
    params:
        chrom = "{chromosome}"
    shell:
        """
        python windowed-insertion-excess.py {input.bam} <(grep ^{params.chrom} {input.bed}) > {output} 2> {log}
        """


rule combine_chromosomes:
    input:
        expand("insertion_excess/{{sample}}-{chromosome}.ie",
               chromosome=CHROMOSOMES)
    output:
        "success-{sample}.txt"
    shell:
        """
        echo "yieehaa" > {output}
        """
