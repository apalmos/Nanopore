cat<<'EOT'>> 4.5_dnadiff.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 4.5_dnadiff

nucmer \
../hg38.fa \
../consensus/output.fa

EOT

samtools bam2fq mapped.sorted.bam | seqtk seq -A > output.fa

cat<<'EOT'>> 4.5_dnadiff.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 4.5_dnadiff

dnadiff -p ref \
../hg38.fa \
../consensus/output.fa

EOT

cat<<'EOT'>> 4.5_nucdiff.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 4.5_nucdiff

nucdiff ../hg38.fa ../consensus/output.fa ./ test

EOT
