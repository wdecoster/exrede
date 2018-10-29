remapped_common-snps.vcf: CG heterozygous SNPs in 50kb upstream region
NCBI remap to GRCh38


picard AddOrReplaceReadGroups RGSM=sample RGLB=ONT RGPL=PROM RGPU=PROM I=d5945_minimap2-large-c9-locus.bam O=d5945-rg.bam
picard AddOrReplaceReadGroups RGSM=sample RGLB=ONT RGPL=PROM RGPU=PROM I=d6843_minimap2-large-c9-locus.bam O=d6843-rg.bam


whatshap phase -o phased-d6843.vcf remapped_common-snps.vcf d6843-rg.bam
whatshap phase -o phased-d5945.vcf remapped_common-snps.vcf d5945-rg.bam
whatshap haplotag -o tagged-d5945.bam phased-d5945.vcf.gz d5945-rg.bam
whatshap haplotag -o tagged-d6843.bam phased-d6843.vcf.gz d6843-rg.bam


samtools view upstream-region-tagged-d5945.bam | cut -f 12- | tr "\t" "\n"  | grep  "^HP:"  | cut -d ':' -f 3 | sort | uniq | while read S; do samtools view -h upstream-region-tagged-d5945.bam |  awk -v tag="HP:i:$S" '($0 ~ /^@/ || index($0,tag)>0)' > ${S}.sam ; done
samtools view 1.sam -o d5945-H1.bam
samtools view 2.sam -o d5945-H2.bam
samtools view upstream-region-tagged-d6843.bam | cut -f 12- | tr "\t" "\n"  | grep  "^HP:"  | cut -d ':' -f 3 | sort | uniq | while read S; do samtools view -h upstream-region-tagged-d6843.bam |  awk -v tag="HP:i:$S" '($0 ~ /^@/ || index($0,tag)>0)' > ${S}.sam ; done
samtools view 1.sam -o d6843-H1.bam
samtools view 2.sam -o d6843-H2.bam

ls d*-H?.bam | parallel 'samtools view {} | cut -f1 | {.}.readIds'

ls d5945-H?.readIds | parallel 'zgrep -F -f {} d5945.summary.txt.gz | cut -f1 > {.}.fofn' &
ls d6843-H?.readIds | parallel 'zgrep -F -f {} d6843.summary.txt.gz | cut -f1 > {.}.fofn'

ls d5945-H?.fofn | parallel 'find /complgen3/ori/promethion/20180302_d5945_wouter/ | grep -F -f {} > {.}.paths' &
ls d6843-H?.fofn | parallel 'find /complgen3/ori/promethion/20180207_wouter_d6843_2/ | grep -F -f {} > {.}.paths' &

for i in $(ls *.paths | sed 's/.paths//')
 do
     mkdir ${i}_fast5
 done

cat d5945-H1.paths | parallel 'cp {} d5945-H1_fast5'
cat d5945-H2.paths | parallel 'cp {} d5945-H2_fast5'
cat d6843-H1.paths | parallel 'cp {} d6843-H1_fast5'
cat d6843-H2.paths | parallel 'cp {} d6843-H2_fast5'

samtools merge all_reads.bam *.bam
samtools index all_reads.bam
samtools fastq all_reads.bam > all_reads.fastq

nanopolish index -d d5945-H1_fast5/ -d d5945-H2_fast5/ -d d6843-H1_fast5/ -d d6843-H2_fast5/ all_reads.fastq
nanopolish call-methylation --reads all_reads.fastq --bam all_reads.bam --genome ~/databases/Homo_sapiens/GRCh38_recommended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -t16 > methylation.txt
