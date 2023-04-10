library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

# Load study table
load("../data/processed/study_table.RData")

# Read in metadata
metadata <- read.csv("../data/processed/precalculated.tsv", sep = "\t")

# Load aa_table
load("../data/processed/nuc_table.RData")


vGenes <- LymphoSeq2::geneFreq(nuc_table, locus = "V", family = TRUE)
vGenes <- vGenes %>% merge(select(metadata, c("sample_name", "group_label")), 
                 by.x = "repertoire_id", by.y = "sample_name")
vGenes$group_label <- factor(vGenes$group_label,      # Reordering group factor levels
                         levels = c("No NACT", "Short Interval", "Long Interval"))
multicolors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Set1")))(28)
p1 <- ggplot(vGenes, aes(x = repertoire_id, y = gene_frequency, fill = gene_name)) +
  geom_bar(stat = "identity") +
  facet_grid(~group_label, scales = "free_x", space = "free_x") +
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0)) + 
  guides(fill = ggplot2::guide_legend(ncol = 2)) +
  scale_fill_manual(values = multicolors) + 
  labs(y = "Frequency (%)", x = "", fill = "") +
  theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
p1
