#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --partition=hive7d
#SBATCH --ntasks=64 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=64 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=256G ## memory to allocate for your job in MB
#SBATCH --time=3-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
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


for R1 in *.MergeBamAlignment.merged.bam; do
 BASENAME="$(basename "${R1%.MergeBamAlignment.merged.bam}")"
 OUT="./$BASENAME"

 echo "${R1} ${BASENAME}"

 gatk --java-options "-Xmx256G -XX:ParallelGCThreads=64" MarkDuplicates \
			--INPUT ${OUT}.MergeBamAlignment.merged.bam \
			--OUTPUT ${OUT}.MarkDuplicates.dedupped.bam \
			--CREATE_INDEX true \
			--VALIDATION_STRINGENCY SILENT \
			--METRICS_FILE ${OUT}.MarkDuplicates.metrics \
                        --TMP_DIR ./temp

done
