#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o de_novo_minimap

directory=$(pwd)
file=$1

minimap2 -X -t 10 -x ava-ont ${directory}/QC/${file} ${directory}/QC/${file} > ${directory}/assembly/AONT.paf

echo 'done de-novo'

#awk '$5=="-"' ${directory}/assembly/AONTself.paf | awk '{print $1}'| sort|uniq > ${directory}/assembly/invertedrepeat.list

miniasm -f ${directory}/QC/${file} ${directory}/assembly/AONT.paf > ${directory}/assembly/AONT.gfa

grep '^S' ${directory}/assembly/AONT.gfa | awk '{print ">"$2"\n"$3}' > ${directory}/assembly/AONT_miniasm.fasta

minimap2 -x asm20 -a --MD -t 8 ${directory}/human_g1k_v37.fasta ${directory}/assembly/AONT_miniasm.fasta > ${directory}/assembly/AONT.sam

#awk '{if($1==$6){print}}' ${directory}/assembly/AONT_miniasm.paf > ${directory}/assembly/AONT_miniasm_self.paf

module load samtools

######output SAM file to BAM file
samtools view -S -b ${directory}/assembly/AONT.sam > ${directory}/assembly/mapped.bam

echo 'finished making bam'

#####sort the BAM file
samtools sort ${directory}/assembly/mapped.bam -o ${directory}/assembly/mapped.sorted.bam

echo 'finished sorting bam'

#####index the BAM file
samtools index ${directory}/assembly/mapped.sorted.bam > ${directory}/assembly/mapped.sorted.bam.bai

echo 'finished indexing bam'

#####see what the BAM file is like
samtools flagstat ${directory}/assembly/mapped.sorted.bam

echo 'run NanoStat for read info'

NanoStat --bam ${directory}/assembly/mapped.sorted.bam > ${directory}/assembly/read_report.txt

echo 'make a vcf'

############ Get a VCF of the alignments
sniffles \
--input ${directory}/assembly/mapped.sorted.bam \
--vcf ${directory}/assembly/variants2.vcf \
--reference ../Resources/human_g1k_v37.fasta
