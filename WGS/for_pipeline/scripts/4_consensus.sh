cat<<'EOT'>> 4_consensus.sh
#!/bin/bash
#SBATCH -n 20
#SBATCH --mem-per-cpu=19G
#SBATCH -t 48:00:00
#SBATCH -o 4_consensus

directory=`cat directory.txt`
echo 'start alignment'

minimap2 \
-a \
--MD \
-x asm10 \
-t 8 \
${directory}/hg38.fa \
${directory}/assembly/dbg.cns.fa \
> ${directory}/consensus/aln.sam

module load samtools

######output SAM file to BAM file
samtools view -S -b ${directory}/consensus/aln.sam > ${directory}/consensus/mapped.bam

echo 'finished making bam'

#####sort the BAM file
samtools sort ${directory}/consensus/mapped.bam -o ${directory}/consensus/mapped.sorted.bam

echo 'finished sorting bam'

#####index the BAM file
samtools index ${directory}/consensus/mapped.sorted.bam

echo 'finished indexing bam'

#####see what the BAM file is like
samtools flagstat ${directory}/consensus/mapped.sorted.bam

echo 'run NanoStat for read info'

NanoStat --bam ${directory}/consensus/mapped.sorted.bam > ${directory}/consensus/read_report.txt

echo 'make a vcf'

sniffles \
--input ${directory}/alignment/mapped.sorted.bam \
--vcf ${directory}/alignment/variants.vcf \
--reference ${directory}/hg38.fa

EOT

sbatch -p cpu ${directory}/consensus/4.5_dnadiff.sh
