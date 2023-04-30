library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(forcats)

# Read in metadata
metadata <- read.csv("data/processed/precalculated.tsv", sep = "\t")
metadata$survivor <- "Long Term"
metadata$survivor[metadata$survival_months < 60] <- "Short Term"
metadata$survivor <- factor(metadata$survivor, levels=c("Short Term", "Long Term"))
sample_names <- metadata$sample_name

# Load aa table
load(file = "data/processed/aa_table.RData")
load(file = "data/processed/study_table.RData")
study_table <- merge(study_table, metadata %>% select(sample_name, group_label), by.x = "repertoire_id", by.y = "sample_name")
study_table <- study_table[is.finite(study_table$junction_aa_length), ]
write.csv(study_table, file = "data/processed/TCRs.csv")

# Declaring sample sets for later use
NoNACT.IDs <- study_table %>% filter(group_label == "No NACT") %>% 
  pull(repertoire_id) %>% unique()
ShortNACT.IDs <- study_table %>% filter(group_label == "Short Interval") %>% 
  pull(repertoire_id) %>% unique()
LongNACT.IDs <- study_table %>% filter(group_label == "Long Interval") %>% 
  pull(repertoire_id) %>% unique()
ShortSurvivor.IDs <- metadata %>% filter(survivor == "Short Term") %>% 
  pull(sample_name)
LongSurvivor.IDs <- metadata %>% filter(survivor == "Long Term") %>% 
  pull(sample_name)

# For any given sample ID, this function returns any cancer associated clusters
getCancerClusters <- function(name) {
  clusters <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  if ("cluster" %in% names(clusters)) {
    return(clusters %>% pull(cluster))
  } else {
    return(vector())
  }
}

# For any given sample ID and set of clusters, this function returns a table
# with the number of TCRs per cluster
getNumTCRsFromCluster <- function(name, clusters) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% 
    merge(TCRs, by = "junction_aa") %>% filter(cluster %in% clusters)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  return(table)
}

# This function calculates a CDF for a given sample
getCDF <- function(name) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% merge(TCRs)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  return(ecdf(table$n))
}

# This function returns the TCRs in each cluster for a given sample
getTCRsPerCluster <- function(name) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% merge(TCRs)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  return(table$n)
}

# Average number of clusters
countClusters <- function(name) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  return(length(unique(TCRs$cluster)))
}
metadata$clusters <- lapply(sample_names, countClusters) %>% unlist

# Average number of cancer clusters
countCancerClusters <- function(name) {
  TCRs <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  return(length(unique(TCRs$cluster)))
}
metadata$cancerclusters <- lapply(sample_names, countCancerClusters) %>% unlist

# Average number of TCRs per cluster
metadata$TCRsPerCluster <- lapply(sample_names, getTCRsPerCluster) %>% lapply(median) %>% unlist

# Proportion Cancer TCRs
# This function calculates the proportion of the repertoire that is cancer associated
calculatePropCancer <- function(name) {
  clusters <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  cancerClusters <- vector()
  if ("cluster" %in% names(clusters)) {
    cancerClusters <- clusters %>% pull(cluster)
  }
  
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% 
    merge(TCRs, by = "junction_aa") %>% filter(cluster %in% cancerClusters)
  countCancerTCRs <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    pull(n) %>% sum()
  
  total <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% pull(duplicate_count) %>% sum
  return(countCancerTCRs / total)
}
metadata$propCancer <- lapply(sample_names, calculatePropCancer) %>% unlist

# Number Cancer TCRs
calculateCountCancer <- function(name) {
  clusters <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  cancerClusters <- vector()
  if ("cluster" %in% names(clusters)) {
    cancerClusters <- clusters %>% pull(cluster)
  }
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% 
    merge(TCRs, by = "junction_aa") %>% filter(cluster %in% cancerClusters)
  countCancerTCRs <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    pull(n) %>% sum()
  return(countCancerTCRs)
}
metadata$CancerTCRs <- lapply(sample_names, calculateCountCancer) %>% unlist

# Compare size of clusters
calculateClusterSizeCancer <- function(name) {
  clusters <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  cancerClusters <- vector()
  if ("cluster" %in% names(clusters)) {
    cancerClusters <- clusters %>% pull(cluster)
  }
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% 
    merge(TCRs, by = "junction_aa") %>% filter(cluster %in% cancerClusters)
  countCancerTCRs <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    pull(n) %>% median()
  return(countCancerTCRs)
}
metadata$cancerTCRsPerCluster <- lapply(sample_names, calculateClusterSizeCancer) %>% unlist

save(metadata, file = "data/processed/metadata.RData")
