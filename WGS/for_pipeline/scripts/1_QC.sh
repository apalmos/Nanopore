cat<<'EOT'>> 1_QC.sh
#!/bin/bash
#SBATCH -n 20
#SBATCH --mem-per-cpu=19G
#SBATCH -t 48:00:00
#SBATCH -o 1_QC

file=$1
directory=`cat directory.txt`

gunzip -c ${directory}/fastq_pass/${file} | NanoFilt \
-l 500 \
--headcrop 10 \
> ${directory}/QC/output/${file}_nanofilt_trimmed.fastq

cat ${directory}/QC/output/${file}_nanofilt_trimmed.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${directory}/QC/output/${file}_nanofilt_trimmed.fa

python ./QAlign/qalign/main.py convert \
--input_fasta ${directory}/QC/output/${file}_nanofilt_trimmed.fa \
--outdir ./QC/qcheck_out/${file} \
--qlevel 2 \
--rc 1 \
--kmerpath ./QAlign/qalign/kmermap_models/r9.4_6mer_nucleotide_template_model.txt

cat ${directory}/QC/qalign/*/* > nanofilt.fastq.fa

EOT


cat */* > nanofilt.fastq.fa

python ./QAlign/qalign/main.py convert \
--input_fasta  \
--outdir ./QC/qcheck_out/${file} \
--qlevel 2 \
--rc 1 \
--kmerpath ./QAlign/qalign/kmermap_models/r9.4_6mer_nucleotide_template_model.txt


ls ./testing/ > test.txt
sed "/\.fa$/d" test.txt > temp && mv temp test.txt

while IFS=$'\t' read -r -a myArray

do

  trimmedfile="$(echo -e "${myArray[1]}" | tr -d '[:space:]')"


  sh test.sh $trimmedfile ${myArray[0]}


done < test.txt

EOT

cat<<'EOT'>> test.sh
#!/bin/bash
#SBATCH -n 20
#SBATCH --mem-per-cpu=19G
#SBATCH -t 48:00:00

file=$1

python ./QAlign/qalign/main.py convert \
--input_fasta hg38.fa \
--outdir ./QC/ \
--qlevel 2 \
--rc 0 \
--kmerpath ./QAlign/qalign/kmermap_models/r9.4_6mer_nucleotide_template_model.txt

EOT


cat */res* > forward.fa
