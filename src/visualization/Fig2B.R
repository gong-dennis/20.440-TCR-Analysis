library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(ggsignif)
library(ggseqlogo)

# Read in metadata
load(file = "data/processed/metadata.RData")
metadata$group_label <- factor(metadata$group_label, levels = c("No NACT", "Short Interval", "Long Interval"))

# Create jitterplot with mean overlay and statistical significance
ggplot(metadata, aes(x = group_label, y = cancerTCRsPerCluster, color = group_label)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(x = "", y = "Median Cancer TCRs Per Cluster") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face="bold")) +
  geom_signif(comparisons = list(c("Short Interval", "No NACT")),
              map_signif_level = TRUE, color = "black")
