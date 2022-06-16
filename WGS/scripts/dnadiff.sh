#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 4.5_dnadiff

directory=`cat directory.txt`
file=$1

nucmer \
../Resources/human_g1k_v37.fasta \
../consensus/${file}

dnadiff -p ref \
../Resources/human_g1k_v37.fasta \
../consensus/${file}

nucdiff \
../Resources/human_g1k_v37.fasta \
../consensus/${file}
