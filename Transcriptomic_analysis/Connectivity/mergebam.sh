#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --partition=hive7d
#SBATCH --ntasks=64 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=64 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=256G ## memory to allocate for your job in MB
#SBATCH --time=6-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
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
ALIGNED_DIR="/lustre1/home/mass/eskalon/hiba/sp_align/bam"


for UNMAPPED in ./*.FastqToSam.unmapped.rg.bam; do
    R1=$(basename "$UNMAPPED" ".FastqToSam.unmapped.rg.bam")
    ALIGNED="${ALIGNED_DIR}/${R1}Aligned.sortedByCoord.out.bam"
    OUT="./${R1}.MergeBamAlignment.merged.bam"

    if [[ ! -f "$ALIGNED" ]]; then
        echo "WARNING: aligned BAM not found for $R1: $ALIGNED"
        continue
    fi
    echo "Unaligned: $UNMAPPED"
    echo "Aligned:   $ALIGNED"
    echo "Output:    $OUT"

    gatk --java-options "-Xmx256G -XX:ParallelGCThreads=64" MergeBamAlignment \
                	--REFERENCE_SEQUENCE ${REF} \
                	--UNMAPPED_BAM ${UNMAPPED} \
                	--ALIGNED_BAM ${ALIGNED}\
                	--OUTPUT ${OUT} \
                	--INCLUDE_SECONDARY_ALIGNMENTS false \
                	--VALIDATION_STRINGENCY SILENT --TMP_DIR ./temp
done
