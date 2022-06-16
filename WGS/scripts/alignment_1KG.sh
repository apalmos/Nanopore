#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o alignment_1KG

directory=$(pwd)
file=$1

############ Align to reference using the de-novo genome
minimap2 \
-a \
--MD \
-x map-ont \
-t 10 \
../Resources/human_g1k_v37.fasta \
${directory}/QC/${file} \
> ${directory}/alignment/aln.sam

############ Samtools to get a cleaned file
module load samtools

######output SAM file to BAM file
samtools view -S -b ${directory}/alignment/aln.sam > ${directory}/alignment/mapped.bam

echo 'finished making bam'

#####sort the BAM file
samtools sort ${directory}/alignment/mapped.bam -o ${directory}/alignment/mapped.sorted.bam

echo 'finished sorting bam'

#####index the BAM file
samtools index ${directory}/alignment/mapped.sorted.bam > ${directory}/alignment/mapped.sorted.bam.bai

echo 'finished indexing bam'

#####see what the BAM file is like
samtools flagstat ${directory}/alignment/mapped.sorted.bam

echo 'run NanoStat for read info'

NanoStat --bam ${directory}/alignment/mapped.sorted.bam > read_report.txt

NanoStat --bam mapped.sorted.bam > read_report.txt


echo 'make a vcf'

############ Get a VCF of the alignments
sniffles \
--input ${directory}/assembly/wtdbg2//mapped.sorted.bam \
--vcf ${directory}/assembly/wtdbg2//variants.vcf \
--reference ../Resources/human_g1k_v37.fasta \
