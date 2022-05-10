cat<<'EOT'>> start_QC.sh
#!/bin/bash
#SBATCH -n 20
#SBATCH --mem-per-cpu=19G
#SBATCH -t 48:00:00

date=`cat ../date.txt`

ls ../../fastq_pass > all.txt
mkdir output

while IFS=$'\t' read -r -a myArray

do

  trimmedfile="$(echo -e "${myArray[1]}" | tr -d '[:space:]')"


  sh QC.sh $trimmedfile ${myArray[0]}


done < all.txt
EOT

cat<<'EOT'>> QC.sh
#!/bin/bash
#SBATCH -n 20
#SBATCH --mem-per-cpu=19G
#SBATCH -t 48:00:00

file=$1
date=`cat ../date.txt`

##porechop -i ../../fastq_pass/${file} \
##-o ./output/${file} \
##--discard_middle \
##--threads 8

gunzip -c ../../fastq_pass/${file} | NanoFilt \
-l 500 \
--headcrop 10 \
> ./output/${file}_nanofilt_trimmed.fastq
gzip ./output/${file}_nanofilt_trimmed.fastq

EOT
