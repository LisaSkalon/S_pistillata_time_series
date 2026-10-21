#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=8 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=8 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=64G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=%J.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=%J.output ## stdout file name(The %J will print job ID number)
#SBATCH --mail-type=NONE ## Send your job status via e-mail: Valid type values are NONE, BEGIN, END, FAIL, REQUEUE, ALL
########### Job information #############
echo "================================"
echo "Start at `date`"
echo "Job id is $SLURM_JOBID"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NTASKS processors."
echo "================================"
#########################################

######## Load required modules ##########
#. /etc/profile.d/modules.sh # Required line for modules environment to work
#module load openmpi/1.8.4 python/2.7 # Load modules that are required by your program
#########################################
source /lustre1/home/mass/eskalon/miniconda/bin/activate gatk2
### Below you can enter your program job command ###


REF="GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.fna"
OUT="GVCFall"


# Extract variant quality scores and make diagnostics plots 
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" SelectVariants \
     -R ${REF} -V  cohort_genotypes.vcf.gz -O ${OUT}.SNPs.vcf.gz \
      -select-type SNP 1> ${OUT}_SNPs.vcf.gz.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" SelectVariants \
     -R ${REF} -V  cohort_genotypes.vcf.gz -O ${OUT}.INDELs.vcf.gz \
     --select-type-to-include INDEL 1> ${OUT}_INDELs.vcf.gz.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantsToTable\
    --variant ${OUT}.SNPs.vcf.gz \
    --output ${OUT}_SNPs.table \
    -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum \
    -F ReadPosRankSum 1> ${OUT}_SNPs.table.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantsToTable \
     --variant ${OUT}.INDELs.vcf.gz \
     --output ${OUT}_INDELs.table -F CHROM -F POS -F QUAL -F QD -F DP \
     -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
     1> ${OUT}_INDELs.table.log 2>&1\

# plot diagnostics plots with these results - script diagnostics.R 
