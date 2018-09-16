### dependencies
# samtools
# bedtools
# jvarkit bamslicebed
# python pysam
# python numpy
# script windowed-insertion-excess.py

set -eo pipefail

usage()
{
cat << EOF
usage: $0 options
This script will execute a pipeline to detect repeat expansions in long read sequencing data
OPTIONS:
   -h | --help             Show this message
   -b | --bam		   bam file to work on
   -f | --fasta		   fasta genome to which the bam is aligned
   -w | --window	   size of window to use
EOF
}

for arg in "$@"; do
  shift
  case "$arg" in
   "--help")    set -- "$@" "-h" ;;
   "--bam")    set -- "$@" "-b"   ;;
   "--fasta")      set -- "$@" "-f"   ;;
   "--window")    set -- "$@" "-w"   ;;
   *)           set -- "$@" "$arg" ;;
  esac
done
# getopts assigns the arguments to variables
while getopts "hb:f:w:" OPTION ; do
  case $OPTION in
   b) BAM=$OPTARG   ;;  # appends an array rather than set a single variable
   f) FASTA=$OPTARG    ;;
   w) WINDOW=$OPTARG  ;;
   h) usage && exit 1 ;;
  esac
done

FAI=${FASTA}.fai
BED=windows-${WINDOW}.bed
BAMSLICED=sliced_${BAM}


samtools faidx $FASTA
bedtools makewindows -g <(cut -f1-2 $FAI) -w $WINDOW | awk -F $'\t' '{printf("%s\t%s\t%d\n",$1,$2,int($3)-1);}' > $BED
samtools index $BAM

for CHROM in $(cut -f1 $FAI)
 do
 bamslicebed -B <(grep "^$CHROM" $BED) $BAM | samtools sort -o $BAMSLICED
 python windowed-insertion-excess.py $BAMSLICED $BED  > insertion_excess_${CHROM}.txt
 done
