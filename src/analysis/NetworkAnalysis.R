library(GGally)
library(network)
library(sna)
library(ggplot2)
library(expm)
library(tidyverse)

# Load in resources
load(file = "data/processed/metadata.RData")
allSampleCluster <- read.csv("data/processed/AllSampleCluster.csv")
load(file = "data/processed/aa_table.RData")

NoNACT.IDs <- metadata %>% filter(group_label == "No NACT") %>% 
  pull(sample_name) %>% unique()
ShortNACT.IDs <- metadata %>% filter(group_label == "Short Interval") %>% 
  pull(sample_name) %>% unique()
LongNACT.IDs <- metadata %>% filter(group_label == "Long Interval") %>% 
  pull(sample_name) %>% unique()
ShortSurvivor.IDs <- metadata %>% filter(survivor == "Short Term") %>% 
  pull(sample_name)
LongSurvivor.IDs <- metadata %>% filter(survivor == "Long Term") %>% 
  pull(sample_name)

# Find the stringent cancer associated TCRs
test <- merge(allSampleCluster, aa_table, by = "junction_aa")

stringentCancerTCRs <- test %>% group_by(cluster) %>% summarize(stringentCancer = sum(stringentCancer == "Yes"),
                                                                stringentNormal = sum(stringentNormal == "Yes"),
                                                                count = sum(duplicate_count), 
                                                                repertoires = n_distinct(repertoire_id),
                                                                unique_aa = n_distinct(junction_aa)) %>%
  filter(stringentCancer > stringentNormal) %>% arrange(-stringentCancer) %>%
  filter(stringentNormal == 0) %>% filter(repertoires > 1)

allSpecificClusters <- test %>% filter(cluster %in% stringentCancerTCRs$cluster) %>% pull(junction_aa)
specificClusters <- stringentCancerTCRs$cluster

# Generate Network Map
networkMap <- test %>% filter(cluster %in% specificClusters)
networkMap <- networkMap %>% merge(select(metadata, c(sample_name, group_label)), by.x = "repertoire_id", by.y = "sample_name")

##### Based on samples #####

# Create dataframe with rows samples and columns clusters
df <- networkMap %>% group_by(cluster, repertoire_id) %>% summarize(duplicate_count=sum(duplicate_count))
cluster_mat <- pivot_wider(data = df, names_from = "cluster", values_from = "duplicate_count") %>% as.data.frame()
cluster_mat[is.na(cluster_mat)] <- 0
rownames(cluster_mat) <- cluster_mat$repertoire_id; cluster_mat$repertoire_id <- NULL

# Create relational matrix
n_samples <- dim(cluster_mat)[1]
my_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    shared_clusters <- sum(sqrt(cluster_mat[i,] * cluster_mat[j,]))
    if (shared_clusters > 0) {
      my_matrix[i, j] <- shared_clusters
    }
  }
}
colnames(my_matrix) <- rownames(cluster_mat)
rownames(my_matrix) <- rownames(cluster_mat)

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
groupNames <- rownames(my_matrix)
groupNames[rownames(my_matrix) %in% NoNACT.IDs] <- "No NACT"
groupNames[rownames(my_matrix) %in% ShortNACT.IDs] <- "Short Interval"
groupNames[rownames(my_matrix) %in% LongNACT.IDs] <- "Long Interval"

survivalNames <- rownames(my_matrix)
survivalNames[rownames(my_matrix) %in% ShortSurvivor.IDs] <- "Short Term"
survivalNames[rownames(my_matrix) %in% LongSurvivor.IDs] <- "Long Term"


net %v% "group" = groupNames
net %v% "survival" = survivalNames

networkMap$colorS <- "#FF9F80"
networkMap$colorS[survivalNames == "Long Term"] <- "#C1E1C5"
net %v% "colorS" = networkMap$colorS
  
networkMap$color <- "#E69F00"
networkMap$color[groupNames == "Short Interval"] <- "#56B4E9"
networkMap$color[groupNames == "Long Interval"] <- "#009E73"
net %v% "color" = networkMap$color

# Plot network
ggnet2(net, color = "color", size = "degree", alpha = 0.4, 
       edge.color = "grey", edge.size = "weight", mode = 'kamadakawai') +
  guides(size = FALSE)

networkCoords <- ggnet2(net, color = "color", size = "degree", alpha = 0.4, 
       edge.color = "grey", edge.size = "weight", mode = 'kamadakawai')$data

degree(net)[groupNames == "Short Interval"] %>% mean(na.rm = TRUE) /2
rowSums(cluster_mat)[groupNames == "Short Interval"] %>% mean(na.rm = TRUE)

degree(net)[groupNames == "Long Interval"] %>% mean(na.rm = TRUE) /2
rowSums(cluster_mat)[groupNames == "Long Interval"] %>% mean(na.rm = TRUE)

degree(net)[groupNames == "No NACT"] %>% mean(na.rm = TRUE) /2
rowSums(cluster_mat)[groupNames == "No NACT"] %>% mean(na.rm = TRUE)

# Plot network
ggnet2(net, color = "colorS", size = "degree", alpha = 0.4, 
       edge.color = "grey", edge.size = "weight", mode = 'kamadakawai') +
  guides(size = FALSE)

networkCoords <- ggnet2(net, color = "colorS", size = "degree", alpha = 0.4, 
                        edge.color = "grey", edge.size = "weight", mode = 'kamadakawai')$data

degree(net)[survivalNames == "Long Term"] %>% mean(na.rm = TRUE) /2
rowSums(cluster_mat)[survivalNames == "Long Term"] %>% mean(na.rm = TRUE)

degree(net)[survivalNames == "Short Term"] %>% mean(na.rm = TRUE) /2
rowSums(cluster_mat)[survivalNames == "Short Term"] %>% mean(na.rm = TRUE)
