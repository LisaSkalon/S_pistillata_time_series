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

# --- 2nd-pass filtering (Genotype Quality) ---

DP_MIN=8.00

# 2nd-pass filtering
gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantFiltration \
 --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED.vcf" \
 --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf" \
            --genotype-filter-name "DP_filter" \
 --genotype-filter-expression "DP < $DP_MIN" \
            1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=8" VariantFiltration \
 --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED.vcf" \
 --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf" \
            --genotype-filter-name "DP_filter" \
 --genotype-filter-expression "DP < $DP_MIN" \
            1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1

# Set filtered sites to no call
gatk SelectVariants --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
    1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1
gatk SelectVariants --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
    1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1


# Extract and plot DP info for each sample (combine it all together and plot; we dont care about each sample as we just want to check that filtering worked)
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table" -F CHROM -F POS -GF GT -GF DP \
    1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table.log" 2>&1
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table" -F CHROM -F POS -GF GT -GF DP \
    1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table.log" 2>&1

# Corrected Loop for SNPs: Using >> to append all samples to one file for the diagnostic plot
SNP_DP_FILE="${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt"
printf "DP\n" > "$SNP_DP_FILE"
while read i;
do
        cut -f $i,$((i+1)) "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk 'NR>1 && $1 != "./." && $1 != ".|." {print $2}' >> "$SNP_DP_FILE"
done < <(head -n 1 "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')

# Corrected Loop for INDELs: Fixed filename and used >> to append
INDEL_DP_FILE="${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt"
printf "DP\n" > "$INDEL_DP_FILE"
while read i;
do
        cut -f $i,$((i+1)) "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk 'NR>1 && $1 != "./." && $1 != ".|." {print $2}' >> "$INDEL_DP_FILE"
done < <(head -n 1 "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.DP.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')

## ## Check number of non-missing genotype DP entries after 2nd-pass filtering
for F in *.GT.DP.txt; do awk 'NR>1' $F; done | wc -l
# 84696172
awk 'NR>1' "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 7866554
awk 'NR>1' "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 716062



# Make diagnostic plots with plot.dp.scores.sh script (Rscript plot_DP_scores.R "${OUT}.DP.afterFiltering")

# 3rd-pass filtering

gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" -F CHROM -F POS -GF GT -GF AD -GF DP \
    1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table.log" 2>&1


# Each sample has 3 columns (GT, AD, DP); get the GT column index for each sample first, extract idx - idx+2, filter using these two columns.
while read i;
do
        GT=$(cut -f $i "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | head -n 1)
        # Corrected: Generating unique filenames per sample to avoid overwriting
        cut -f $i-$((i+2)) "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '$1 != "./." && $1 != ".|." { if(NR==1) {print $2} else {split($2,a,","); for (j in a) {print a[j]/$3} } }' > $GT.AD.txt
done < <(head -n 1 "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')

# Make diagnostic plots with plot.dp.scores.sh script (Rscript plot_AD_scores.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD" and Rscript plot_AD_scores.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.XaxisLimits" )

# Generate allelic depth plots using sum(AD) as denominator instead of DP
while read i;
do
        GT=$(cut -f $i "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | head -n 1)
        # Corrected: Generating unique filenames per sample to avoid overwriting
        cut -f $i-$((i+2)) "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '$1 != "./." && $1 != ".|." { if(NR==1) {print $2} else {split($2,a,","); SUM=0; for (j in a) {SUM=SUM+a[j]}; if(SUM>0){for (j in a) {print a[j]/SUM}} } }' > $GT.AD.sumSDdenominator.txt
done < <(head -n 1 "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" | awk '{for (i=1; i<=NF; ++i) {if($i~".GT$") {print i}}}')

# Make diagnostic plots with plot.dp.scores.sh script (Rscript plot_AD_scores.sumSDdenominator.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.sumSDdenominator" andRscript plot_AD_scores.sumSDdenominator.R "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.sumSDdenominator.XaxisLimits")

# VCF to Table (use for phylogeny)
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf"   --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.table"   -F CHROM -F POS -GF GT \
    1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.table.log" 2>&1
gatk VariantsToTable --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.table" -F CHROM -F POS -GF GT \
    1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.table.log" 2>&1


