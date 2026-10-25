library(igraph)
library(ggraph)
library(ggplot2)
library(stringr)
library(dplyr)


#save.image(file='kinship5.RData')
#load('kinship.RData')
#### all clones/no relatives ####
setwd("/home/gospozha/haifa/hiba/sp_align/snp/pruned")
# Read KING table
king <- read.csv("idb.kin0", stringsAsFactors=FALSE, sep = "\t")
#king <- read.csv("idb_0.7.kin0", stringsAsFactors=FALSE, sep = "\t") # all
#king <- read.csv("./new/idb_w_rel.kin0", stringsAsFactors=FALSE, sep = "\t") # no_relatives
#king <- read.csv("idb.kin0.no_known_clones.tsv", stringsAsFactors=FALSE, sep = "\t") # no_relatives

meta <- read.csv("Metadata.csv", stringsAsFactors=FALSE)

# Make IDs comparable between KING and metadata
king$ID1 <- sub("^\\./", "", king$ID1)
king$ID2 <- sub("^\\./", "", king$ID2)
meta$id <- as.character(meta$id)

# Define thresholds
clone_cutoff <- 0.354      # definite clones
firstdeg_cutoff <- 0.177   # first-degree relatives

king$rel_type <- ifelse(king$KINSHIP > clone_cutoff, "Clone",
                        ifelse(king$KINSHIP > firstdeg_cutoff, "1st-degree", NA))

# Keep clone + first-degree relationships
edges <- subset(king, !is.na(rel_type), select=c(ID1, ID2, KINSHIP, rel_type))

# Build graph
g <- graph_from_data_frame(edges, directed=FALSE)

# Add metadata to graph vertices
vertex_meta <- data.frame(name=V(g)$name, stringsAsFactors=FALSE) %>%
  left_join(meta, by=c("name"="id"))

# Check whether any graph samples failed to match Metadata.csv
vertex_meta %>% filter(is.na(depth))


# add metadata
V(g)$month <- vertex_meta$month
V(g)$date <- vertex_meta$data
V(g)$depth <- factor(vertex_meta$depth, levels=c(5, 10, 25, 40))
V(g)$species <- vertex_meta$species
V(g)$label <- vertex_meta$name

# plot
kinship <- ggraph(g, layout="fr") +
  geom_edge_link(aes(color=rel_type, width=KINSHIP), alpha=0.6) +
  geom_node_point(aes(color=depth, shape=species), size=4) +
  geom_node_text(aes(label=label), repel=TRUE, size=3) +
  scale_edge_color_manual(values=c("Clone"="#6F2DBD", "1st-degree"="#F4A261")) +
  scale_color_manual(values=c("5"="#E64B35", "10"="#4DBBD5", "25"="#00A087", "40"="#3C5488")) +
  scale_edge_width_continuous(name="Kinship", range=c(0.3, 2)) +
  theme_void() +
  labs(edge_color="Relationship", color="Depth (m)", shape="Species")

kinship

ggsave("kinship.jpg", kinship, width=5, height=5, dpi=300)

#saveRDS(kinship, "fig_kinship.rds")

edges_meta <- edges %>%
  left_join(meta, by=c("ID1"="id")) %>%
  rename(month1=month, date1=data, depth1=depth, species1=species) %>%
  left_join(meta, by=c("ID2"="id")) %>%
  rename(month2=month, date2=data, depth2=depth, species2=species)

View(edges_meta)

write.csv(edges_meta, "clones_table.csv")