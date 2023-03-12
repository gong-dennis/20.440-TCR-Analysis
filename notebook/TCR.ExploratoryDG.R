library(tidyverse)
library(LymphoSeq2)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

study_table <- LymphoSeq2::readImmunoSeq("../data/ChemoProjTCRs/")
study_statistics <- read.table("../data/SampleOverview.tsv", sep = '\t')

summary_table <- LymphoSeq2::clonality(study_table)
summary_table