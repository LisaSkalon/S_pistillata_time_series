#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=8 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=8 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=64G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=GVCFall.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=GVCFall.output ## stdout file name(The %J will print job ID number)
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


# 1st-pass filtering - values are chosed according to the diagnostics plots

# SNP Filter Thresholds
SNP_QD_MIN=2.0
SNP_MQ_MIN=40.0
SNP_FS_MAX=60.0
SNP_SOR_MAX=4.0

# Indel Filter Thresholds
INDEL_QD_MIN=2.0
INDEL_MQ_MIN=40.0
INDEL_FS_MAX=200.0
INDEL_SOR_MAX=10.0

# General Quality
QUAL=30.0


gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" \
        VariantFiltration --reference ${REF} --variant ${OUT}.SNPs.vcf.gz  \
        --output ${OUT}_SNPs_VarScores_filter.qd.vcf.gz \
        -filter "QUAL < $QUAL"  --filter-name "VarScores_filter_QUAL"  \
        -filter  "QD < $SNP_QD_MIN"  --filter-name "VarScores_filter_QD" \
        -filter  "MQ < $SNP_MQ_MIN" --filter-name "VarScores_filter_MQ"  \
        -filter "FS > $SNP_FS_MAX"  --filter-name "VarScores_filter_FS"  \
        -filter  "SOR > $SNP_SOR_MAX" --filter-name "VarScores_filter_SOR" \
        1> ${OUT}_SNPs_VarScores_filter.qd.vcf.gz.log 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" \
        VariantFiltration --reference ${REF} --variant ${OUT}.INDELs.vcf.gz \
        --output ${OUT}_INDELs_VarScores_filter.qd.vcf.gz \
        -filter "QUAL < $QUAL"  --filter-name "VarScores_filter_QUAL"  \
	-filter  "QD < $INDEL_QD_MIN"  --filter-name "VarScores_filter_QD" \
	-filter  "MQ < $INDEL_MQ_MIN" --filter-name "VarScores_filter_MQ"  \
	-filter "FS > $INDEL_FS_MAX"  --filter-name "VarScores_filter_FS"  \
	-filter  "SOR > $INDEL_SOR_MAX" --filter-name "VarScores_filter_SOR" \
 	1> ${OUT}_INDELs_VarScores_filter.qd.vcf.gz.log 2>&1

# Check the number of PASSED after the first filtering
zcat "${OUT}.SNPs.vcf.gz" | grep -v '^#' | wc -l
# 1433670
zcat "${OUT}_SNPs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l
# 1175072
zcat "${OUT}.INDELs.vcf.gz" | grep -v '^#' | wc -l
# 96273
zcat "${OUT}_INDELs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l
# 95398


# Extract only variants that PASSED filtering
zcat "${OUT}_SNPs_VarScores_filter.qd.vcf.gz" | grep -E '^#|PASS' > "${OUT}_SNPs_VarScores_filterPASSED.vcf"
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" IndexFeatureFile --input "${OUT}_SNPs_VarScores_filterPASSED.vcf" 1> "${OUT}_SNPs_VarScores_filterPASSED.vcf.log" 2>&1
zcat "${OUT}_INDELs_VarScores_filter.qd.vcf.gz" | grep -E '^#|PASS' > "${OUT}_INDELs_VarScores_filterPASSED.vcf"
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" IndexFeatureFile --input "${OUT}_INDELs_VarScores_filterPASSED.vcf" 1> "${OUT}_INDELs_VarScores_filterPASSED.vcf.log" 2>&1



# Extract and plot DP info for each sample (from all samples before previous filtering; shows us overall variant coverage, independent of other factors)
gatk VariantsToTable --reference "${REF}" --variant cohort_genotypes.vcf.gz --output "${OUT}.DP.table" -F CHROM -F POS -GF GT -GF DP 1> "${OUT}.DP.table.log" 2>&1

# Each sample has 2 columns (GT, DP); get the GT column index for each sample first, extract idx and idx+1, filter using these two columns.
while read i;
do
        GT=$(cut -f $i "${OUT}.DP.table" | head -n 1)
        cut -f $i,$((i+1)) "${OUT}.DP.table" | awk '$1 != "./." && $1 != ".|." {print $2}' > $GT.DP.txt
done < <(head -n 1 "${OUT}.DP.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')

# Make diagnostic plots and choose the apropriate DP value with plot.dp.scores.sh script (Rscript plot_DP_scores.R "${OUT}.DP")
