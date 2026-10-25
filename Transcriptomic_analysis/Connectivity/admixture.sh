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
#source /lustre1/home/mass/eskalon/miniconda/bin/activate bcftools
### Below you can enter your program job command ###

# convert scaffolds to integer for admixture
cp VCF_pruned_30.bim VCF_pruned_30.bim.original

awk 'BEGIN{OFS="\t"}
{
    if (!($1 in chr)) {
        n++
        chr[$1]=n
    }
    $1=chr[$1]
    print
}' VCF_pruned_30.bim.original > VCF_pruned_30.bim

# run admixture
for K in 1 2 3 4 5 6 7 8 9 10; do
    for SEED in 41 42 43 44 45; do
        ./dist/admixture_linux-1.3.0/admixture \
            --cv=10 \
            -j8 \
            -s "$SEED" \
            VCF_pruned_30.bed "$K" \
            2>&1 | tee "admixture.K${K}.seed${SEED}.log"
    done
done

# compare CV across K and seeds
echo -e "K\tSeed\tCV_error" > admixture_CV.tsv

for f in admixture.K*.seed*.log; do
    K=$(echo "$f" | sed -E 's/.*K([0-9]+)\.seed.*/\1/')
    SEED=$(echo "$f" | sed -E 's/.*seed([0-9]+)\.log/\1/')
    CV=$(grep "CV error" "$f" | awk '{print $NF}')

    echo -e "${K}\t${SEED}\t${CV}" >> admixture_CV.tsv
done

{ head -n 1 admixture_CV.tsv; tail -n +2 admixture_CV.tsv | sort -k1,1n -k2,2n; } \
    > admixture_CV_sorted.tsv
