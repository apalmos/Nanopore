cat<<'EOT'>> alignment_pipeline.sh
#!/bin/bash
#SBATCH -n 20
#SBATCH --mem-per-cpu=19G
#SBATCH -t 48:00:00

echo 'start alignment'

minimap2 \
-a \
--MD \
-x map-ont \
-R '@RG\tID:PAI07137\tSM:C_dil' \
-t 8 \
/scratch/groups/sgdp_nanopore/Resources/hg38.fa \
../error_correction/medaka_output/consensus.fasta \
> aln.sam

module load apps/samtools

######output SAM file to BAM file
samtools view -S -b aln.sam > mapped.bam

echo 'finished making bam'

#####sort the BAM file
samtools sort mapped.bam -o mapped.sorted.bam

echo 'finished sorting bam'

#####index the BAM file
samtools index mapped.sorted.bam

echo 'finished indexing bam'

#####see what the BAM file is like
samtools flagstat mapped.sorted.bam

echo 'run NanoStat for read info'

NanoStat --bam mapped.sorted.bam > read_report.txt

echo 'make a vcf'

sniffles \
--input mapped.sorted.bam \
--vcf variants.vcf \
--reference /scratch/groups/sgdp_nanopore/Resources/hg38.fa


EOT

/scratch/groups/sgdp_nanopore/software/bedops/bedops/bin/convert2bed \
--input=vcf \
--do-not-sort \
--deletions \
< variants.vcf > deletions.bed

#/scratch/groups/sgdp_nanopore/software/bedops/bedops/bin/bedextract \
#chr9 \
#insertions.bed > chr9.bed

###/scratch/groups/sgdp_nanopore/software/bedops/bedops/bin/sort-bed chr9.bed > test

EOT
