gatk HaplotypeCaller \
   -R genome.fa \
   -I sample.sort.dedup.realign.bam \
   -O sample.raw_variants.vcf

samtools bam2fq input.bam | seqtk seq -A > output.fa

./work/d1/3e5a7bb4f36f358ac42156fca41753/variants.vcf

samtools faidx genome.fa
samtools dict genome.fa > genome.dict

cat<<'EOT'>> gatk_indel.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00

gatk SelectVariants \
-R  /scratch/prj/sgdp_nanopore/Resources/human_g1k_v37.fasta \
-V variants.vcf \
-O raw_indels.vcf \
--select-type-to-include INDEL

EOT

gatk SelectVariants \
-R  /scratch/prj/sgdp_nanopore/Resources/human_g1k_v37.fasta \
-V variants.vcf \
-O raw_snps.vcf \
--select-type-to-include SNP

gatk VariantsToTable \
  -V variants.vcf \
  -F CHROM -F POS -F TYPE -GF AD \
  -O output.table

vcftools --vcf variants.vcf --chr 1 --from-bp 1000000 --to-bp 2000000

vcftools --vcf variants.vcf --freq --out output

vcftools --vcf variants.vcf --depth -c > depth_summary.txt

gatk Funcotator \
     --variant variants.vcf \
     --reference /scratch/prj/sgdp_nanopore/Resources/human_g1k_v37.fasta \
     --ref-version hg19 \
     --data-sources-path /scratch/prj/sgdp_nanopore/Resources/funcotator_dataSources.v1.7.20200521s \
     --output variants.funcotated.vcf \
     --output-file-format VCF


gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download




gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
    -V variants.vcf \
    --trust-all-polymorphic \
    --tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode INDEL \
    --max-gaussians 4 \
    --resource:mills,known=false,training=true,truth=true,prior=12: "/scratch/prj/sgdp_nanopore/Resources/Mills_and_1000G_gold_standard.indels.hg38.vcf" \
    --resource:axiomPoly,known=false,training=true,truth=false,prior=10: "/scratch/prj/sgdp_nanopore/Resources/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf" \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2: "/scratch/prj/sgdp_nanopore/Resources/Homo_sapiens_assembly38.dbsnp138.vcf" \
    -O cohort_indels.recal \
    --tranches-file cohort_indels.tranches




gatk ValidateVariants -V /scratch/prj/sgdp_nanopore/Resources/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf



   gatk MakeSitesOnlyVcf \
           -I variants.vcf \
           -O cohort_sitesonly.vcf.gz


samtools addreplacerg -r "@RG\tID:S1" -o test.bam sortedbam.bam


samtools view -H test.bam | grep '@RG'

tabix -p vcf variants.vcf


gatk BaseRecalibrator \
-R /scratch/prj/sgdp_nanopore/Resources/hg38.fa \
-I test.bam \
-known-sites variants.vcf \
-O sample_1_recal_data.table


gatk Funcotator \
     --variant variants.vcf \
     --reference /scratch/prj/sgdp_nanopore/Resources/hg38.fa \
     --ref-version hg38 \
     --data-sources-path /scratch/prj/sgdp_nanopore/Resources/funcotator_dataSources.v1.7.20200521g \
     --output variants.funcotated.maf \
     --output-file-format MAF


     gatk SelectVariants \
     -R  /scratch/prj/sgdp_nanopore/Resources/hg38.fa \
     -V variants.vcf \
     -O raw_snps.vcf \
     --select-type-to-include SNP



bcftools query -f '%CHROM %POS %ID %REF %ALT %QUAL %FILTER %AF\n' variants.vcf | head -3


bcftools query -f '%CHROM %POS %ID %REF %ALT %QUAL %AF\n' variants.vcf > test.vcf

snpEff download -c ./snpEff/snpEff.config -v hg38

snpEff eff -c ./snpEff/snpEff.config hg38 test.vcf > test2.vcf




#NanoPlot
NanoPlot --summary /scratch/prj/sgdp_nanopore/20_01_22/Alish_sample1/20220118_1558_1F_PAI08255_70f005aa/sequencing_summary_PAI08255_02964990.txt --loglength -o summary-plots-log-transformed

NanoPlot -t 12 --color yellow --bam sortedbam.bam --downsample 10000 -o bamplots_downsampled






#mapped reads
samtools view -b -F 4 file.bam > mapped.bam

#unmapped reads
samtools view -b -f 4 file.bam > unmapped.bam





bcftools view -f PASS variants.vcf > test.vcf



cat<<'EOT'>> samtools.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00

bcftools mpileup -f /scratch/prj/sgdp_nanopore/Resources/hg38.fa sortedbam.bam | bcftools call -mv -Ob -o calls.bcf

EOT




samtools mpileup -g -f /scratch/prj/sgdp_nanopore/Resources/hg38.fa sortedbam.bam > myraw.bcf

bcftools view -bvcg my-raw.bcf > my-var.bcf

bcftools view my.var.bcf |
vcfutils.pl varFilter - > my.var-final.vcf

cat<<'EOT'>> phase.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00

whatshap phase -o phased.vcf --ignore-read-groups --reference=/scratch/prj/sgdp_nanopore/Resources/hg38.fa variants.vcf sortedbam.bam

EOT

srun -p cpu --pty /bin/bash

cat<<'EOT'>> duet.sh
#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00

duet \
sortedbam.bam \
/scratch/prj/sgdp_nanopore/Resources/hg38.fa \
./duet

EOT

conda install -c bioconda clair3 python=3.6.10 -y
conda create --name coloc --clone seurat4
conda remove --name seurat4 --all

conda install -c conda-forge/label/gcc7 r-rmysql

cat<<'EOT'>> RStudio.sh
#!/bin/bash -l
#SBATCH --job-name=test-rstudio
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --signal=USR2
#SBATCH --cpus-per-task=1

# get unused socket per https://unix.stackexchange.com/a/132524
export PASSWORD=$(openssl rand -base64 15)
readonly IPADDRESS=$(hostname -I | tr ' ' '\n' | grep '10.211.4.')
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:
   ssh -J $USER@bastion.er.kcl.ac.uk -NL 8787:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk
   and point your web browser to http://localhost:8787
2. Login to RStudio Server using the following credentials:
   user: ${USER}
   password: ${PASSWORD}
When done using the RStudio Server, terminate the job by:
1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:
      scancel -f ${SLURM_JOB_ID}
END

export RSTUDIOBASEDIR=$HOME/rstudiosbase
export RSTUDIOTMPDIR=$RSTUDIOBASEDIR/tmp
export RSTUDIORUNDIR=$RSTUDIOBASEDIR/run
export RSTUDIOVARLIBDIR=RSTUDIOBASEDIR/var-lib
mkdir -p $RSTUDIOTMPDIR $RSTUDIORUNDIR $RSTUDIOVARLIBDIR

singularity exec --bind /scratch:/scratch --bind $RSTUDIOTMPDIR:/tmp --bind $RSTUDIORUNDIR:/run --bind $RSTUDIOVARLIBDIR:/var/lib/rstudio-server /software/containers/singularity/rstudio/rstudio_4.1.1.sif rserver --server-user ${USER} --www-port ${PORT} --auth-none=0 --auth-pam-helper-path=pam-helper
printf 'RStudio Server exited' 1>&2

EOT

./configure --prefix=/users/k1463257/R/x86_64-pc-linux-gnu-library/4.1/libxml2-2.9.14
