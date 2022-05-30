cat<<'EOT'>> 2_de_novo.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 2_de_novo

directory=`cat directory.txt`
file=$1

minimap2 -X -t 10 -x ava-ont ${directory}/QC/nanofilt.fastq.gz ${directory}/QC/nanofilt.fastq.gz > ${directory}/assembly/AONT.paf

echo 'done de-novo'

awk '{if($1==$6){print}}' ${directory}/assembly/AONT.paf > ${directory}/assembly/AONTself.paf

awk '$5=="-"' ${directory}/assembly/AONTself.paf | awk '{print $1}'| sort|uniq > invertedrepeat.list

miniasm -f ${directory}/QC/nanofilt.fa ${directory}/assembly/AONT.paf > ${directory}/assembly/AONT.gfa

grep '^S' ${directory}/assembly/AONT.gfa | awk '{print ">"$2"\n"$3}' > ${directory}/assembly/AONT_miniasm.fasta

minimap2 -x asm10 ${directory}/assembly/AONT_miniasm.fasta ${directory}/assembly/AONT_miniasm.fasta -t 1 -X > ${directory}/assembly/AONT_miniasm.paf

awk '{if($1==$6){print}}' ${directory}/assembly/AONT_miniasm.paf > ${directory}/assembly/AONT_miniasm_self.paf

echo 'done assembly'

echo 'done making fasta'

EOT


cat<<'EOT'>> 2_de_novo.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 2_de_novo

directory=`cat directory.txt`

./wtdbg2/wtdbg2 -x ont -g 300m -i \
${directory}/QC/nanofilt.fa \
-t 16 -fo ${directory}/assembly/dbg

./wtdbg2/wtpoa-cns -t 16 -i ${directory}/assembly/dbg.ctg.lay.gz -fo ${directory}/assembly/dbg.raw.fa

module load samtools

minimap2 -t16 -ax map-ont -r2k \
${directory}/assembly/dbg.raw.fa \
${directory}/QC/nanofilt.fa | samtools sort -@4 > ${directory}/assembly/dbg.bam

samtools view -F0x900 ${directory}/assembly/dbg.bam | ./wtdbg2/wtpoa-cns -t 16 -d ${directory}/assembly/dbg.raw.fa -i - -fo ${directory}/assembly/dbg.cns.fa

bwa index ${directory}/assembly/dbg.cns.fa

EOT

cat<<'EOT'>> 2_de_novo.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 2_de_novo

directory=`cat directory.txt`

module load samtools

minimap2 -t16 -ax map-ont -r2k \
${directory}/assembly/dbg.raw.fa \
${directory}/QC/nanofilt.fa | samtools sort -@4 > ${directory}/assembly/dbg.bam

samtools view -F0x900 ${directory}/assembly/dbg.bam | ./wtdbg2/wtpoa-cns -t 16 -d ${directory}/assembly/dbg.raw.fa -i - -fo ${directory}/assembly/dbg.cns.fa

bwa index ${directory}/assembly/dbg.cns.fa

EOT




cat<<'EOT'>> miniasm.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o miniasm

directory=`cat directory.txt`

miniasm \
-c 2 \
-h 1100 \
-d 7000 \
-f \
${directory}/QC/nanofilt.fa \
${directory}/assembly/minimap.paf > ${directory}/assembly/assemble.gfa

EOT


awk 'NF' forward.fa > test && mv test forward.fa
awk 'NF' reverse.fa > test && mv test reverse.fa
