library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(forcats)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

# Read in metadata
metadata <- read.csv("../data/processed/precalculated.tsv", sep = "\t")

# Load aa table
load(file = "../data/processed/aa_table.RData")
load(file = "../data/processed/study_table.RData")
study_table <- merge(study_table, metadata %>% select(sample_name, group_label), by.x = "repertoire_id", by.y = "sample_name")
study_table <- study_table[is.finite(study_table$junction_aa_length), ]

# Declaring sample sets for later use
NoNACT.IDs <- study_table %>% filter(group_label == "No NACT") %>% 
  pull(repertoire_id) %>% unique()
ShortNACT.IDs <- study_table %>% filter(group_label == "Short Interval") %>% 
  pull(repertoire_id) %>% unique()
LongNACT.IDs <- study_table %>% filter(group_label == "Long Interval") %>% 
  pull(repertoire_id) %>% unique()

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

write.csv(study_table, file = "TCRs.csv")

##### Analysis from 4/23 #####
setwd("~/Documents/GitHub/20.440-TCR-Analysis/")


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


##### Creation of stacked bar chart #####

# Calculate prevalence of cancer associated TCRs and plot


caTCRPrevalenceTable <- function(name) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% merge(TCRs)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  TCRs <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  table$cancer <- table$cluster %in% TCRs$cluster
  table <- head(table, 20)
  table$prop <- table$n / sum(table$n)
  table$sample_name <- name
  return(arrange(table, -prop))
}

hi <- lapply(sample_names, caTCRPrevalenceTable) %>% do.call(what = rbind)
hi <- hi %>% merge(select(metadata, c(sample_name, group_label)), by = "sample_name")
hi$cancer[hi$cancer == TRUE] <- "Cancer Associated"
hi$cancer[hi$cancer == FALSE] <- "Non-Cancer Associated"

my_colors <- c("#9AC9E3", "#EFC3E6")
ggplot(hi, aes(x = sample_name, y = prop,
               fill = fct_reorder(factor(cancer), prop),
               color = cancer)) +
  geom_col(position = "fill", colour="black", size = 0.2) +
  facet_grid(~ group_label, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "", y = "Proportion of 20 Top Clusters", fill = "") +
  theme_classic() + ylim(0,1) +
  theme(legend.position = "top") +
  scale_fill_manual(values = my_colors)


##### Clustering all samples together #####


allSampleCluster <- read.csv("data/processed/AllSampleCluster.csv")

test <- merge(allSampleCluster, aa_table, by = "junction_aa")

test_plot_stringent <- test %>% group_by(cluster, stringentCancer, stringentNormal) %>% summarize(Distinct_IDs = n_distinct(repertoire_id), n = n())


test_plot$stringent <- "Broadly Expressed"
test_plot$stringent[test_plot$cluster %in% stringentCancerTCRs$cluster] <- "Cancer Specific"

stringentCancerTCRs <- test %>% group_by(cluster) %>% summarize(stringentCancer = sum(stringentCancer == "Yes"),
                                         stringentNormal = sum(stringentNormal == "Yes")) %>%
  filter(stringentCancer > stringentNormal) %>% arrange(-stringentCancer)

test_plot <- test %>% group_by(cluster, cancer) %>% summarize(Distinct_IDs = n_distinct(repertoire_id), n = n()) %>% arrange(-Distinct_IDs)

# Cancer Specific
cluster6063 <- test %>% filter(cluster == 6063) %>% pull(junction_aa)
cluster3454 <- test %>% filter(cluster == 3454) %>% pull(junction_aa)

allSpecificClusters <- test %>% filter(cluster %in% stringentCancerTCRs$cluster) %>% pull(junction_aa)
specificClusters <- stringentCancerTCRs$cluster

# Broadly expressed
cluster5380 <- test %>% filter(cluster == 5380) %>% pull(junction_aa)



