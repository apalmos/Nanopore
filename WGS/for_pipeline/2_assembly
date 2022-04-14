cat<<'EOT'>> de_novo.sh
#!/bin/bash
#SBATCH -n 4
#SBATCH --mem-per-cpu=9G
#SBATCH -t 48:00:00

minimap2 \
-x ava-ont -t 16 \
/mnt/lustre/groups/sgdp_nanopore/20_01_22/analysis/QC/nanofilt.gz \
/mnt/lustre/groups/sgdp_nanopore/20_01_22/analysis/QC/nanofilt.gz \
> minimap.paf

echo 'done de-novo'

miniasm -f \
../QC/nanofilt.fastq.gz \
minimap.paf > assemble.gfa

echo 'done assembly'

awk '/^S/{print ">"$2"\n"$3}' assemble.gfa > miniasm.fasta

echo 'done making fasta'

EOT



minimap2 -x ava-ont -t 16 \
../QC/nanofilt.fastq.gz \
../QC/nanofilt.fastq.gz \
> minimap.paf


| gzip -1 > ./minimap.paf.gz


############ Get a report of the data
assembly-stats ./miniasm.fasta > report.txt

############ Run DNADiff
dnadiff -p dnadiff \
/scratch/groups/sgdp_nanopore/Resources/hg38.fa \
./miniasm.fasta

############ Run MUMMMER
mummer \
/scratch/groups/sgdp_nanopore/Resources/hg38.fa \
./miniasm.fasta

############ Align to reference using the de-novo genome
/scratch/groups/sgdp_nanopore/software/minimap2/minimap2 \
-a \
--MD \
-x map-ont \
-t 8 \
/scratch/groups/sgdp_nanopore/Resources/hg38.fa \
./miniasm.fasta
> aln.sam

############ Samtools to get a cleaned file
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

############ Get a VCF of the alignments
/scratch/groups/sgdp_nanopore/software/Sniffles-master/bin/sniffles-core-1.0.12/sniffles \
-m mapped.sorted.bam \
-v variants.vcf

EOT
