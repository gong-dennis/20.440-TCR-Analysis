library(GGally)
library(network)
library(sna)
library(ggplot2)
library(expm)
library(tidyverse)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

metadata <- read.csv("../data/processed/precalculated.tsv", sep = "\t")
allSampleCluster <- read.csv("../data/processed/AllSampleCluster.csv")
load(file = "../data/processed/aa_table.RData")

NoNACT.IDs <- metadata %>% filter(group_label == "No NACT") %>% 
  pull(sample_name) %>% unique()
ShortNACT.IDs <- metadata %>% filter(group_label == "Short Interval") %>% 
  pull(sample_name) %>% unique()
LongNACT.IDs <- metadata %>% filter(group_label == "Long Interval") %>% 
  pull(sample_name) %>% unique()

test <- merge(allSampleCluster, aa_table, by = "junction_aa")

stringentCancerTCRs <- test %>% group_by(cluster) %>% summarize(stringentCancer = sum(stringentCancer == "Yes"),
                                                                stringentNormal = sum(stringentNormal == "Yes")) %>%
  filter(stringentCancer > stringentNormal) %>% arrange(-stringentCancer)
allSpecificClusters <- test %>% filter(cluster %in% stringentCancerTCRs$cluster) %>% pull(junction_aa)
specificClusters <- stringentCancerTCRs$cluster

networkMap <- test %>% filter(cluster %in% specificClusters)
networkMap <- networkMap %>% merge(select(metadata, c(sample_name, group_label)), by.x = "repertoire_id", by.y = "sample_name")

##### Based on samples #####

df <- networkMap %>% group_by(cluster, repertoire_id) %>% summarize(duplicate_count=sum(duplicate_count))
df_wide <- pivot_wider(data = df, names_from = "cluster", values_from = "duplicate_count")
df_wide <- as.data.frame(df_wide)
df_wide[is.na(df_wide)] <- 0
rownames(df_wide) <- df_wide$repertoire_id
df_wide$repertoire_id <- NULL
cluster_mat <- df_wide

my_matrix <- sqrt(as.matrix(cluster_mat) %*% t(as.matrix(cluster_mat)))
#my_matrix[my_matrix > 0] <- log2(my_matrix[my_matrix > 0])
my_matrix[my_matrix > 0] <- 3 * my_matrix[my_matrix > 0] / max(my_matrix)

my_matrix[lower.tri(my_matrix)] <- my_matrix[upper.tri(my_matrix)]
diag(my_matrix) <- 0

# Initialize the network
net <- as.network(x = my_matrix, # the network object
                  directed = FALSE, # specify whether the network is directed
                  loops = FALSE, # do we allow self ties (should not allow them)
                  matrix.type = "adjacency",# the type of input
                  names.eval = "weight",
                  ignore.eval = FALSE)

# Network properties
network.vertex.names(net) <- rownames(my_matrix)
groupNames <- metadata$group_label
groupNames[rownames(my_matrix) %in% NoNACT.IDs] <- "No NACT"
groupNames[rownames(my_matrix) %in% ShortNACT.IDs] <- "Short Interval"
groupNames[rownames(my_matrix) %in% LongNACT.IDs] <- "Long Interval"

net %v% "group" = groupNames
networkMap$color <- "#E69F00"
networkMap$color[groupNames == "Short Interval"] <- "#56B4E9"
networkMap$color[groupNames == "Long Interval"] <- "#009E73"
net %v% "color" = networkMap$color

# Plot network
ggnet2(net, color = "color", size = "degree", alpha = 0.4, 
       edge.color = "grey", edge.size = "weight", mode = 'princoord') +
  guides(size = FALSE)

networkCoords <- ggnet2(net, color = "color", size = "degree", alpha = 0.4, 
       edge.color = "grey", edge.size = "weight", mode = 'kamadakawai')$data
