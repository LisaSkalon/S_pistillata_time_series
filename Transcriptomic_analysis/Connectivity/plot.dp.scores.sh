#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=1 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=1 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=12G ## memory to allocate for your job in MB
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
source /lustre1/home/mass/eskalon/miniconda/bin/activate r-env
### Below you can enter your program job command ###
DP_MIN=20.000
OUT="GVCFall"

#Rscript plot_DP_scores.R "${OUT}.DP" $DP_MIN 100 1> "${OUT}.DP.log" 2>&1
Rscript plot_DP_scores.R "${OUT}.DP.afterFiltering" $DP_MIN 100 1> "${OUT}.DP.afterFiltering.log" 2>&1

Rscript plot_AD_scores.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD" 0.0 1.0 1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.log" 2>&1
Rscript plot_AD_scores.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.XaxisLimits" 0.1 0.9 1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.XaxisLimits.log" 2>&1

Rscript plot_AD_scores.sumSDdenominator.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.sumSDdenominator" 0.0 1.0 1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.sumSDdenominator.log" 2>&1
Rscript plot_AD_scores.sumSDdenominator.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.sumSDdenominator.XaxisLimits" 0.1 0.9 1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.sumSDdenominator.XaxisLimits.log" 2>&1

