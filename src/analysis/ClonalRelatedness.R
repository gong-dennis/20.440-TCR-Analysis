library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

# Load study table
load("../data/processed/study_table.RData")

# Read in metadata
metadata <- read.csv("../data/processed/precalculated.tsv", sep = "\t")

# Clonal relatedness

clonalRelatedness <- LymphoSeq2::clonalRelatedness(study_table = study_table, edit_distance = 10)

test <- metadata %>% select(c("sample_name", "group_label")) %>% 
  merge(aa_table, by.x = "sample_name", by.y = "repertoire_id")
colnames(test)[1] <- "repertoire_id"
p1 <- LymphoSeq2::topSeqsPlot(study_table = filter(test, group_label == "No NACT"), top = 10)
p2 <- LymphoSeq2::topSeqsPlot(study_table = filter(test, group_label == "Short Interval"), top = 10)
p3 <- LymphoSeq2::topSeqsPlot(study_table = filter(test, group_label == "Long Interval"), top = 10)

grid.arrange(
  p1, p2, p3,
  ncol = 3
)