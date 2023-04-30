library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(ggsignif)
library(ggseqlogo)

# Load in data
load(file = "data/processed/aa_table.RData")
allSampleCluster <- read.csv("data/processed/AllSampleCluster.csv")

# Match clusters with amino acids
test <- merge(allSampleCluster, aa_table, by = "junction_aa")

# Cancer Associated
cluster6063 <- test %>% filter(cluster == 6063) %>% pull(junction_aa)

# Cancer Specific
cluster3454 <- test %>% filter(cluster == 3454) %>% pull(junction_aa)

# Broadly expressed
cluster5380 <- test %>% filter(cluster == 5380) %>% pull(junction_aa)

# Logo plots
ggseqlogo(cluster6063, method = 'prob' )
ggseqlogo(cluster3454, method = 'prob' )
ggseqlogo(cluster5380, method = 'prob' )
