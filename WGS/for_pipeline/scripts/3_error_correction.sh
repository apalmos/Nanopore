cat<<'EOT'>> 3_error_correction.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 3_error_correction

directory=`cat directory.txt`

minipolish \
-t 16 \
--rounds 2 \
${directory}/QC/nanofilt.fastq.gz \
${directory}/assembly/assemble.gfa > ${directory}/error_correction/minipolished_assembly.gfa

echo 'done polishing'

awk '/^S/{print ">"$2"\n"$3}' ${directory}/error_correction/minipolished_assembly.gfa | fold > ${directory}/error_correction/polished.fasta

echo 'exported fasta file'

EOT

cat<<'EOT'>> medaka.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o medaka

directory=`cat directory.txt`

medaka_consensus \
-t 14 \
-i ${directory}/QC/nanofilt.fastq.gz \
-d ${directory}/error_correction/polished.fasta \
-o ${directory}/consensus/medaka_output \
–m r941_prom_high_g303

echo 'done making consensus sequence'

EOT

sbatch -p cpu ${directory}/error_correction/4_consensus.sh

bwa index ${directory}/error_correction/polished.fasta

cat<<'EOT'>> bwa.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o bwa

directory=`cat directory.txt`

bwa mem \
${directory}/error_correction/polished.fasta \
${directory}/QC/nanofilt.fastq.gz > ${directory}/error_correction/bwa_mapping.sam

module load samtools

samtools view -Sb ${directory}/error_correction/bwa_mapping.sam > ${directory}/error_correction/bwa_mapping.bam
samtools sort –o b${directory}/error_correction/wa_mapping.sorted.bam ${directory}/error_correction/bwa_mapping.bam
samtools index ${directory}/error_correction/bwa_mapping.sorted.bam

EOT
