## RNA-seq SNPs analysis of _S. pistillata_ samples 
Bash scripts are based on the pipeline for SNP analysis of [Coral_Stress_Phenome](https://github.com/hputnam/Coral_Stress_Phenome/tree/main/Genotype_Analysis/Pocillopora_acuta_PacBio_Assembly/RNAseq_short_variant_analysis) project and on the work of [Federica Scucchia](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs), with modifications.

### Programs
GATK4 pipeline was installed on the remote cluster through conda using this [guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4). PLINK2 and bcftools were also installed via conda. rgsam was installed through [github](https://github.com/djhshih/rgsam), and then the Makefile was manually changed to ensure installation to the local home bin folder.

### Prepare BAM files: Convert reads to BAM, add read group info, merge alignments, filter 
[Fastqtobam.sh](fastqtobam.sh) script converts paired fastq files to BAM files sorted by read name, [rgsam.sh](rgsam.sh) adds read group info (RG) to aligned reads, [dictreference.sh](dictreference.sh) indexes and creates a dictionary from the reference, [mergebam.sh](mergebam.sh) merges unaligned BAM files with aligned BAMs and additionally filters the alinged read (e.g. removes secondary alignments). 

### Mark duplicates
[markdup.sh](markdup.sh) marks potential PCR duplicates.

### Split reads that contain Ns in their CIGAR string 
[splitsigar.sh](splitsigar.sh) splits reads with N in CIGAR - used for post-processing of the alignment.

### Call SNPs and indels
[haplotcall.sh](haplotcall.sh) calls variants (SNPs and indels) simultaneously via local de-novo assembly of haplotypes in an active region. 

### Combine *.g.vcf.gz files and call genotypes
[combinegvcf.sh](combinegvcf.sh) combines all gvcf files to a cohort file, [genotypegvcf.sh](genotypegvcf.sh) performs joint genotyping on a cohort. 

### Filter SNPs and indels
[filtvar1.sh](filtvar1.sh) selects SNPs and indels and extracts variant quality scores for making diagnostics plots with [diagnostics.R](diagnostics.R).
[filtvar2.sh](filtvar2.sh)  performs 1st-pass filtering of variants with the parameters chosen according to the diagnostics plots.
[filtvar3.sh](filtvar3.sh)  performs 2nd and 3rd-pass filtering with the parameters chosen according to the diagnostics plots drawn by [plot.dp.scores.sh](plot.dp.scores.sh) script.   

Number of variants passed:
```
# Check the number of PASSED after the first filtering
zcat "${OUT}.SNPs.vcf.gz" | grep -v '^#' | wc -l
# 1477466
zcat "${OUT}_SNPs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l
# 1078283
zcat "${OUT}.INDELs.vcf.gz" | grep -v '^#' | wc -l
# 123034
zcat "${OUT}_INDELs_VarScores_filter.qd.vcf.gz" | grep 'PASS' | wc -l
# 91994
```
```
# Check number of non-missing genotype DP entries after 2nd-pass filtering
for F in *.GT.DP.txt; do awk 'NR>1' $F; done | wc -l
# 26999368
awk 'NR>1' "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 4462395
awk 'NR>1' "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.GT.DP.txt" | wc -l
# 445852
```
