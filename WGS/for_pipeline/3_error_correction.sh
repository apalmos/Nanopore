cat<<'EOT'>> error_correction.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=36G
#SBATCH -t 48:00:00

minimap2 -t 16 \
../miniasm/miniasm.fasta \
../QC/nanofilt.fastq.gz \
 > minimap.racon.paf

 echo 'done minimap alignment to de-novo sequence'

racon \
../QC/nanofilt.fastq.gz
-t 4 \
minimap.racon.paf \
../miniasm/miniasm.fasta \
> miniasm.racon.consensus.fasta

echo 'done racon'

EOT

cat<<'EOT'>> error_correction3.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=36G
#SBATCH -t 48:00:00

minipolish \
../QC/nanofilt/nanofilt.fastq.gz \
../miniasm/assemble.gfa > minipolished_assembly_t16.gfa

echo 'done polishing'

awk '/^S/{print ">"$2"\n"$3}' minipolished_assembly_t16.gfa | fold > filtered_t16.fasta

echo 'exported fasta file'

EOT

conda activate medaka

medaka_consensus -i filtered.fasta \
-d miniasm.racon.consensus.fasta \
-o medaka_output \
–m r941_prom_high_g303

conda deactivate

echo 'done making consensus sequence'

EOT








cat<<'EOT'>> error_correction.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=36G
#SBATCH -t 48:00:00

racon \
-t 4 \
../QC/nanofilt/nanofilt.fastq.gz \
minimap.racon.paf \
../miniasm/miniasm.fasta \
> ./miniasm.racon.consensus.fasta

EOT

cat<<'EOT'>> error_correction2.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=36G
#SBATCH -t 48:00:00

minipolish \
-t 4 \
../QC/nanofilt/nanofilt.fastq.gz \
../miniasm/assemble.gfa > minipolished_assembly.gfa

echo 'done polishing'

awk '/^S/{print ">"$2"\n"$3}' minipolished_assembly.gfa | fold > filtered.fasta

echo 'exported fasta file'

EOT

cat<<'EOT'>> error_correction3.sh
#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=36G
#SBATCH -t 48:00:00

conda activate medaka

medaka_consensus -i filtered.fasta \
-d miniasm.racon.consensus.fasta \
-o medaka_output \
–m r941_prom_high_g303

conda deactivate

echo 'done making consensus sequence'

EOT
