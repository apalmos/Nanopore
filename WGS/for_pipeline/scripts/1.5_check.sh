#!/bin/bash

directory=`cat directory.txt`

cat ${directory}/QC/output/*.gz > ${directory}/nanofilt.fastq.gz

sbatch -p cpu ${directory}/assembly/de_novo.sh
sbatch -p cpu ${directory}/alignment/align.sh




cat<<'EOT'>> 1.5_qalign.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 1.5_qalign

zcat nanofilt.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > file.fa

rm nanofilt.fastq.gz

cd ./QAlign

python qalign/main.py convert \
--input_fasta ../file.fa \
--outdir ../QC/ \
--qlevel 2 \
--rc 0 \
--kmerpath qalign/kmermap_models/r9.4_6mer_nucleotide_template_model.txt

EOT
