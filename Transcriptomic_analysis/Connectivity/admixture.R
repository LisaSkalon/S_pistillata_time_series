library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

setwd("/home/gospozha/haifa/hiba/sp_align/snp/admixture/")



# ADMIXTURE ancestry coefficients
Q <- read.table("VCF_pruned_30.K5.best.Q")

colnames(Q) <- paste0("Cluster", 1:5)

# Sample order used by ADMIXTURE
fam <- read.table("VCF_pruned_30.fam", stringsAsFactors = FALSE)

samples <- data.frame(
  FID = fam$V1,
  IID = fam$V2
)

# remove "./" so IDs match Metadata.csv
samples$id <- sub("^\\./", "", samples$IID)

# Add ancestry proportions
admix <- cbind(samples, Q)

# Metadata
meta <- read.csv("Metadata.csv", stringsAsFactors = FALSE)
meta$id <- as.character(meta$id)

admix <- admix %>%
  left_join(meta, by = "id")

admix %>%
  filter(is.na(depth))

# long format
admix_long <- admix %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "Cluster",
    values_to = "Ancestry"
  )

# sort by species and depth
admix <- admix %>%
  arrange(species, depth, id) %>%
  mutate(order = row_number())

admix_long <- admix %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "Cluster",
    values_to = "Ancestry"
  )

ggplot(admix_long,
       aes(x = order, y = Ancestry, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(~ species + depth,
             scales = "free_x",
             space = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Individuals",
    y = "Ancestry proportion"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.3, "lines")
  )


#### without separating by depth and species apriori ####
# Assign each individual to the cluster with its highest ancestry proportion
Qcols <- paste0("Cluster", 1:5)

admix$dominant_cluster <- Qcols[
  max.col(admix[, Qcols], ties.method = "first")
]

# Maximum ancestry proportion, useful for sorting within cluster
admix$maxQ <- apply(admix[, Qcols], 1, max)

# Label containing biological information
admix$sample_label <- paste0(admix$species, "_", admix$depth, "m")

# Order first by genetic cluster, then species/depth,
# then by strength of assignment
admix <- admix %>%
  arrange(dominant_cluster, species, depth, desc(maxQ)) %>%
  mutate(plot_order = row_number())

# Long format for plotting
admix_long <- admix %>%
  pivot_longer(
    cols = all_of(Qcols),
    names_to = "Cluster",
    values_to = "Ancestry"
  )

ggplot(admix_long,
       aes(x = plot_order,
           y = Ancestry,
           fill = Cluster)) +
  
  geom_col(width = 1) +
  
  facet_grid(
    ~ dominant_cluster,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_x_continuous(
    breaks = admix$plot_order,
    labels = admix$sample_label,
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "Ancestry proportion"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6
    ),
    panel.spacing.x = unit(0.5, "lines")
  )

#### tables and statistics ####

admix %>%
  count(dominant_cluster, depth) %>%
  group_by(dominant_cluster) %>%
  mutate(percent = 100 * n / sum(n)) -> depth.table

admix %>%
  count(dominant_cluster, species) %>%
  group_by(dominant_cluster) %>%
  mutate(percent = 100 * n / sum(n)) -> species.table

# Make sure variables are factors
admix$dominant_cluster <- factor(admix$dominant_cluster)
admix$species <- factor(admix$species)
admix$depth <- factor(admix$depth, levels = c(5, 10, 25, 40))


# species difference at 40 m
admix_40 <- subset(admix, depth == 40)

tab_species_40 <- table(
  admix_40$dominant_cluster,
  admix_40$species
)

tab_species_40

fisher_species_40 <- fisher.test(
  tab_species_40,
  simulate.p.value = TRUE,
  B = 100000
)

fisher_species_40
# 0.8119
# does genetic-cluster membership differ between pistillata and kuehlmanni?

prop.table(tab_species_40, margin = 2)

# depth association of pistillata
admix_pist <- subset(admix, species == "pistillata")

tab_depth <- table(
  admix_pist$dominant_cluster,
  admix_pist$depth
)

tab_depth

fisher_depth <- fisher.test(
  tab_depth,
  simulate.p.value = TRUE,
  B = 100000
)

fisher_depth
# 1e-05
# is cluster associated with a specific depth
prop.table(tab_depth, margin = 2)

# Cramer V effect size
cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  k <- min(nrow(tab) - 1, ncol(tab) - 1)
  sqrt(as.numeric(chi$statistic) / (n * k))
}

cramers_v(tab_depth) 
# 0.615


# summary tables

admix %>%
  count(dominant_cluster, species, depth) %>%
  arrange(dominant_cluster, species, depth) -> summary.table

admix_pist %>%
  count(depth, dominant_cluster) %>%
  group_by(depth) %>%
  mutate(percent = 100 * n / sum(n)) -> summary.table.percent

write.csv(summary.table, "summary.table.admixture.csv", row.names = F)
write.csv(summary.table.percent, "summary.table.percent.admixture.csv", row.names = F)



# if i need additionally filter by SNP distance
bim <- read.table("VCF_pruned_30.bim.original", header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("chr", "snp", "cm", "pos", "a1", "a2")

# Sort SNPs by scaffold and position
bim <- bim[order(bim$chr, bim$pos), ]

# Distance between consecutive SNPs on the same scaffold
bim$distance <- ave(bim$pos, bim$chr, FUN = function(x) c(NA, diff(x)))

# Summary
summary(bim$distance)
sum(bim$distance < 2500, na.rm = TRUE)
sum(bim$distance < 5000, na.rm = TRUE)
sum(bim$distance < 10000, na.rm = TRUE)

# Number of unique scaffolds
length(unique(bim$chr))

ggplot(bim, aes(x = distance)) +
  geom_histogram(bins = 80, fill = "#4575b4") +
  scale_x_log10() +
  labs(x = "Distance to preceding SNP (bp, log10)", y = "Number of SNPs") +
  theme_classic()