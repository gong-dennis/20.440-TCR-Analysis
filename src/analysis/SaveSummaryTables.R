library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(survminer)

##### Set up #####

# Read in TCR sequences
study_table <- LymphoSeq2::readImmunoSeq("data/raw/ChemoProjTCRs/") %>% 
  topSeqs(top = 1000)
save(study_table, file = "data/processed/study_table.RData")

# Read in metadata
metadata <- read.csv("data/processed/precalculated.tsv", sep = "\t")

# Calculate summary statistics
summary_table <- LymphoSeq2::clonality(study_table)
save(summary_table, file = "data/processed/summary_table.RData")

sample_names <- study_table %>% pull(repertoire_id) %>% unique()

# Extract productive sequences
aa_table <- LymphoSeq2::productiveSeq(study_table = study_table, 
                                      aggregate = "junction_aa", prevalence = TRUE)
save(aa_table, file = "data/processed/aa_table.RData")

nuc_table <- LymphoSeq2::productiveSeq(study_table = study_table, aggregate = "junction", 
                                       prevalence = FALSE)
save(nuc_table, file = "data/processed/nuc_table.RData")