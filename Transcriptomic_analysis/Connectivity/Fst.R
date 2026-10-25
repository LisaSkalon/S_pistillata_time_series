library(gdsfmt)
library(SNPRelate)
library(hierfstat)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(gridGraphics)
library(cowplot)
library(StAMPP)
library(adegenet)
library(ggrepel)
library(poppr)
library(patchwork)
library(ape)
library(ggfun)
library(ggtree)
setwd("/home/gospozha/haifa/hiba/sp_align/snp/admixture/")

## Load the pruned SNPs

bed.fn <- "VCF_pruned_30.bed"
fam.fn <- "VCF_pruned_30.fam"
bim.fn <- "VCF_pruned_30.bim.original"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned.gds", cvt.chr="char")
snpgdsSummary("plink_pruned.gds")
## Open the gds file
genofile <- snpgdsOpen("plink_pruned.gds")
###

# Get sample IDs from GDS
sample.id.gds <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Remove ./ prefix for matching
sample.id <- sub("^\\./", "", as.character(sample.id.gds))

# Load your metadata (assumed to be a CSV with columns: sample, site, depth)
meta <- read.csv("Metadata.csv")
meta$id <- as.character(meta$id)
# Samples in GDS but not in metadata
setdiff(sample.id, meta$id)

# Samples in metadata but not in GDS
setdiff(meta$id, sample.id)

# Ensure sample names match
meta <- meta[match(sample.id, meta$id), ]

# Reset row names
rownames(meta) <- NULL

# Verify that the order is now identical
identical(sample.id, meta$id)


#### stampp ####

# Extract genotypes: individuals x SNPs
gds_data <- snpgdsGetGeno(
  genofile,
  snpfirstdim = FALSE,
  with.id = TRUE
)

geno <- gds_data$genotype

# Sample IDs
sample.id <- sub("^\\./", "", as.character(gds_data$sample.id))

# Load full metadata
meta_all <- read.csv("Metadata.csv", stringsAsFactors = FALSE)

meta_all$id <- as.character(meta_all$id)

# Match metadata to genotype order
idx <- match(sample.id, meta_all$id)

meta <- meta_all[idx, , drop = FALSE]
rownames(meta) <- NULL

# Standardize metadata
meta$species <- tolower(trimws(as.character(meta$species)))
meta$depth <- as.character(meta$depth)

# Inspect sample sizes
table(meta$species, meta$depth)

# Select all pistillata and kuehlmanni from 40 m
sel <- which(meta$species == "pistillata" | (meta$species == "kuehlmanni" & meta$depth == "40"))

# Create population labels
pop <- ifelse(meta$species[sel] == "kuehlmanni", "40m K", paste0(meta$depth[sel], "m P"))
table(pop)

# Select individuals and genotypes
g <- geno[sel, , drop = FALSE]
ids <- sample.id[sel]

stopifnot(nrow(g) == length(pop), length(unique(pop)) >= 2, all(table(pop) >= 2))

# SNPRelate uses 3 for missing genotypes; convert to NA if necessary
g[g == 3] <- NA

# Check SNP call rate within each population
min_callrate <- 0.7

callrate <- sapply(unique(pop), function(p) {
  colMeans(!is.na(g[pop == p, , drop = FALSE]))
})

# Retain SNPs with at least 70% genotyping success in each population
keep <- rowSums(callrate >= min_callrate) == ncol(callrate)

# Remove monomorphic SNPs
allele_freq <- colMeans(g / 2, na.rm = TRUE)
keep <- keep & is.finite(allele_freq) & allele_freq > 0 & allele_freq < 1

g <- g[, keep, drop = FALSE]

cat("Individuals:", nrow(g), "\n")
cat("SNPs retained:", ncol(g), "\n")
stopifnot(ncol(g) >= 2)

# Convert genotype dosage (0/1/2) to allele frequencies (0/0.5/1)
freq <- g / 2
colnames(freq) <- paste0("SNP_", which(keep))

# Prepare StAMPP input
genotype <- data.frame(
  Sample = ids, Pop = pop, Ploidy = 2, Format = "freq",
  as.data.frame(freq), check.names = FALSE
)

# Convert to StAMPP format
genotype.st <- stamppConvert(genotype, "r")

# Calculate pairwise FST
set.seed(42)
fst_combined <- stamppFst(genotype.st, nboots = 1000, percent = 95, nclusters = 3)

# Results
fst_combined$Fsts
fst_combined$Pvalues
#fst_combined$Bootstraps

# Save
write.csv(fst_combined$Fsts, "Fst_combined_Fst.csv")
write.csv(fst_combined$Pvalues, "Fst_combined_Pvalues.csv")
saveRDS(fst_combined, "Fst_combined_StAMPP.rds")
#### plot Fst ####
fst_mat <- as.matrix(fst_combined$Fsts)
p_mat <- as.matrix(fst_combined$Pvalues)

# Make matrices symmetric (StAMPP may return only one triangle)
make_symmetric <- function(m) {
  missing <- is.na(m) & !is.na(t(m))
  m[missing] <- t(m)[missing]
  return(m)
}

fst_mat <- make_symmetric(fst_mat)
p_mat <- make_symmetric(p_mat)

# Reorder populations
custom_order <- c("5m P", "10m P", "25m P", "40m P", "40m K")
fst_mat <- fst_mat[custom_order, custom_order]
p_mat <- p_mat[custom_order, custom_order]

# Matrix for coloring: lower triangle only
fst_colors <- fst_mat
fst_colors[upper.tri(fst_colors, diag = TRUE)] <- NA
fst_colors[!is.na(fst_colors) & fst_colors < 0] <- 0

# Labels: FST below diagonal, P values above
n <- nrow(fst_mat)
label_matrix <- matrix("", n, n)

for (i in 1:n) {
  for (j in 1:n) {
    if (i > j && !is.na(fst_mat[i, j])) {
      label_matrix[i, j] <- sprintf("%.3f", fst_mat[i, j])
    } else if (i < j && !is.na(p_mat[i, j])) {
      label_matrix[i, j] <- ifelse(p_mat[i, j] < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", p_mat[i, j])))
    }
  }
}

# Colors
breaks <- seq(0, max(fst_colors, na.rm = TRUE), length.out = 5)
colors <- c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#5a5eaa")
col_fun <- colorRamp2(breaks, colors)

# plot
ht <- Heatmap(
  fst_colors, name = "Fst", col = col_fun, na_col = "white",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left", column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = NA),
  
  heatmap_legend_param = list(
    at = breaks, labels = sprintf("%.2f", breaks),
    title = "Fst", title_gp = gpar(fontsize = 10, fontface = "plain"),
    labels_gp = gpar(fontsize = 9)
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    # FST labels below diagonal, P values above
    label <- label_matrix[i, j]
    if (label != "") {
      text_col <- ifelse(i < j, "grey40", "black")
      grid.text(label, x, y, gp = gpar(fontsize = 9, col = text_col))
    }
    
    # Diagonal line
    if (i == j) {
      grid.lines(x = unit.c(x - 0.5 * width, x + 0.5 * width),
                 y = unit.c(y + 0.5 * height, y - 0.5 * height),
                 gp = gpar(col = "grey80", lwd = 1))
    }
  }
)

draw(ht)

pdf("Fst_combined_pvalues.pdf", width = 6.5, height = 5.5)
draw(ht)
dev.off()

#### PCA ####

# Extract genotypes and population labels
geno.nas <- genotype.st[, -(1:5)]
pop.vector <- factor(genotype.st$Pop, levels = c("5m P", "10m P", "25m P", "40m P", "40m K"))
sample.ids <- as.character(genotype.st$Sample)

# Convert allele frequencies to genind format
convert_to_genind_format <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 1, "11",
                ifelse(x == 0.5, "12",
                       ifelse(x == 0, "22", NA))))
}

geno.char <- as.data.frame(lapply(geno.nas, convert_to_genind_format), stringsAsFactors = FALSE)
rownames(geno.char) <- sample.ids

geno.obj <- df2genind(geno.char, pop = pop.vector, ploidy = 2, NA.char = NA, sep = "")

# calculate PCA
pca <- dudi.pca(tab(geno.obj, NA.method = "mean"), scannf = FALSE, nf = 3)

# Extract individual scores
scores <- as.data.frame(pca$li)
scores$Pop <- pop.vector
scores$SampleID <- sample.ids

# Variance explained
eig <- 100 * pca$eig / sum(pca$eig)

# Calculate centroids
centroids <- scores %>%
  group_by(Pop) %>%
  summarise(Axis1.centroid = mean(Axis1), Axis2.centroid = mean(Axis2), .groups = "drop")

scores <- scores %>% left_join(centroids, by = "Pop")

# plot PCA
pt_to_mm <- function(pt) pt / 2.845

condition_colors <- c(
  "5m P"  = "#D55E00",
  "10m P" = "#0072B2",
  "25m P" = "#009E73",
  "40m P" = "#8172B3",
  "40m K" = "#CC79A7"
)

pca_plot2 <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  geom_segment(aes(xend = Axis1.centroid, yend = Axis2.centroid), alpha = 0.4, linewidth = 0.3) +
  geom_text(data = centroids, aes(x = Axis1.centroid, y = Axis2.centroid, label = Pop),
            size = pt_to_mm(8), show.legend = FALSE, alpha = 0.8) +
  geom_point(size = 1.5) +
  labs(x = paste0("PC1 (", round(eig[1], 1), "%)"),
       y = paste0("PC2 (", round(eig[2], 1), "%)")) +
  scale_color_manual(values = condition_colors) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))

pca_plot2

ggsave("PCA_Stylophora_depth_species.pdf", pca_plot2, width = 6, height = 5)

#### IBS clustering ####

# Identity by State
ibs <- snpgdsIBS(genofile, autosome.only = FALSE, num.thread = 8)

# Get sample IDs
ibs_ids <- sub("^\\./", "", as.character(ibs$sample.id))

# Convert IBS similarity to genetic distance
rownames(ibs$ibs) <- colnames(ibs$ibs) <- ibs_ids
dist_matrix <- as.dist(1 - ibs$ibs)

# Hierarchical clustering (UPGMA)
standard_hclust <- hclust(dist_matrix, method = "average")

# Convert to phylogenetic tree format
ibs_phylo <- as.phylo(standard_hclust)

# Check labels
head(ibs_phylo$tip.label)

meta_all <- read.csv("Metadata.csv", stringsAsFactors = FALSE)
meta_all$id <- as.character(meta_all$id)

# Match metadata to IBS sample IDs
stopifnot(!anyDuplicated(meta_all$id), !anyDuplicated(ibs_ids))
stopifnot(all(ibs_ids %in% meta_all$id))

meta_ibs <- meta_all[match(ibs_ids, meta_all$id), ]
rownames(meta_ibs) <- NULL

# Prepare group labels
meta_ibs$species <- tolower(trimws(meta_ibs$species))
meta_ibs$depth <- as.character(meta_ibs$depth)

meta_ibs$Group <- ifelse(meta_ibs$species == "kuehlmanni", "40m K", paste0(meta_ibs$depth, "m P"))
meta_ibs$Group <- factor(meta_ibs$Group, levels = c("5m P", "10m P", "25m P", "40m P", "40m K"))

# Prepare labels for ggtree
meta_ibs$label <- meta_ibs$id

# Verify matching and sample distribution
stopifnot(identical(ibs_ids, meta_ibs$id), !anyNA(meta_ibs$Group))
table(meta_ibs$Group)
# Remove existing label column
meta_ibs$label <- NULL

# Create label as the first column
meta_ibs <- data.frame(label = as.character(meta_ibs$id),
                       meta_ibs[, setdiff(names(meta_ibs), "id"), drop = FALSE])

# Check
colnames(meta_ibs)
anyDuplicated(colnames(meta_ibs))
setdiff(ibs_phylo$tip.label, meta_ibs$label)
# plot
condition_colors <- c(
  "5m P" = "#D55E00", "10m P" = "#0072B2", "25m P" = "#009E73",
  "40m P" = "#8172B3", "40m K" = "#CC79A7"
)

ibs_plot <- ggtree(ibs_phylo, layout = "rectangular") %<+% meta_ibs +
  geom_tiplab(size = 2, offset = 0.005) +
  geom_tippoint(aes(color = Group, shape = species), size = 2.5) +
  scale_color_manual(values = condition_colors) +
  scale_shape_manual(values = c("pistillata" = 16, "kuehlmanni" = 17)) +
  theme_tree2() +
  labs(color = "Group", shape = "Species")

ibs_plot

ggsave("IBS_Stylophora.pdf", ibs_plot, width = 8, height = 16)
saveRDS(ibs_plot, "IBS_Stylophora.rds")