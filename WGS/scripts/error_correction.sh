#!/bin/bash
#SBATCH -n 40
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 3_error_correction

directory=$(pwd)
file=$1

minipolish \
-t 16 \
--rounds 2 \
${directory}/QC/${file} \
${directory}/assembly/assemble.gfa > ${directory}/error_correction/minipolished_assembly.gfa

echo 'done polishing'

awk '/^S/{print ">"$2"\n"$3}' ${directory}/error_correction/minipolished_assembly.gfa | fold > ${directory}/error_correction/polished.fasta

echo 'exported fasta file'

medaka_consensus \
-t 14 \
-i ${directory}/QC/${file} \
-d ${directory}/error_correction/polished.fasta \
-o ${directory}/consensus/medaka_output \
–m r941_prom_high_g303

echo 'done making consensus sequence'
