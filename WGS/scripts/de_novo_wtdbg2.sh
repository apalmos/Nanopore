#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o de_novo_wtdbg2

file=dummy_original.fastq.gz
directory=`cat directory.txt`

######assemble reads
./wtdbg2/wtdbg2 -x ont -i \
${directory}/QC/${file} \
-t 16 -R -fo ${directory}/assembly/dbg

######create raw consensus
./wtdbg2/wtpoa-cns -t 16 -f -i ${directory}/assembly/dbg.ctg.lay.gz -fo ${directory}/assembly/dbg.raw.fa

######check the raw consensus
NanoStat --fasta ${directory}/assembly/dbg.raw.fa

module load samtools

######align raw reads to raw consensus
minimap2 -t16 -ax map-ont \
${directory}/assembly/dbg.raw.fa \
${directory}/QC/${file} | samtools sort -@4 > ${directory}/assembly/dbg.bam

######output SAM file to BAM file
NanoStat --bam ${directory}/assembly/dbg.bam

samtools view ${directory}/assembly/dbg.bam | ./wtdbg2/wtpoa-cns -t 16 -d ${directory}/assembly/dbg.raw.fa -fo ${directory}/assembly/dbg.cns.fa

NanoStat --fasta ${directory}/assembly/dbg.cns.fa

bwa index ${directory}/assembly/dbg.cns.fa

minimap2 -x asm20 -a --MD -t 8 \
../Resources/human_g1k_v37.fasta \
 ${directory}/assembly/dbg.cns.fa > ${directory}/assembly/AONT.sam

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
--vcf ${directory}/assembly/variants.vcf \
--reference ../Resources/human_g1k_v37.fasta
