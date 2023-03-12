library(tidyverse)
library(LymphoSeq2)
library(ggplot2)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

# LymphoSeq2
study_table <- LymphoSeq2::readImmunoSeq("../data/raw/ChemoProjTCRs/")
summary_table <- LymphoSeq2::clonality(study_table)
summary_table

# From Adaptive
study_statistics <- read.table("../data/analysis/SampleOverview.tsv", 
                               header = TRUE, sep = '\t')

# Visualize histogram

ggplot(study_statistics, aes(x = productive_simpson_clonality)) +
  geom_histogram()
