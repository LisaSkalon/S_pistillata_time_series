#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=8 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=8 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=64G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=plink.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=plink.output ## stdout file name(The %J will print job ID number)
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
source /lustre1/home/mass/eskalon/miniconda/bin/activate bcftools
### Below you can enter your program job command ###


REF="GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.fna"
OUT="GVCFall"

bcftools annotate \
 --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
 -o GVCFall_SNPs_VarScores_filterPASSED_DPfilterNoCall.annotated.vcf \
 --threads 8 \
 GVCFall_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf

plink2 --vcf GVCFall_SNPs_VarScores_filterPASSED_DPfilterNoCall.annotated.vcf \
 --make-bed --chr-set 28 no-xy --allow-extra-chr --double-id --max-alleles 2 \
 --out VCF_annotated

plink2 --bfile VCF_annotated \
 --allow-extra-chr --indep-pairwise 25 5 0.5

plink2 --vcf GVCFall_SNPs_VarScores_filterPASSED_DPfilterNoCall.annotated.vcf \
 --extract plink2.prune.in --allow-extra-chr --make-bed --double-id \
 --out Spist_vcf_pruned
