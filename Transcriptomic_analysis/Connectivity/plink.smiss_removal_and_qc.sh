#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=1 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=1 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=16G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=plink.clones.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=plink.clones.output ## stdout file name(The %J will print job ID number)
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
plink2 --bfile VCF_noclones \
  --allow-extra-chr \
  --remove smiss_to_remove.txt \
  --make-bed \
  --out VCF_nosmiss

#plink2 --bfile VCF_noclones \
#--missing \
#--allow-extra-chr \
#--out missing.noclones

plink2 --bfile VCF_nosmiss \
  --allow-extra-chr \
  --geno 0.3 \
  --maf 0.05 \
  --make-bed \
  --out VCF_nosmiss30

plink2 --bfile VCF_nosmiss30 \
--missing \
--allow-extra-chr \
--out missing.nosmiss30
