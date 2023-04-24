library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(ggsignif)
library(ggseqlogo)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

# Read in metadata
metadata <- read.csv("../data/processed/precalculated.tsv", sep = "\t")

metadata$group_label <- factor(metadata$group_label, levels = c("No NACT", "Short Interval", "Long Interval"))

# Create jitterplot with mean overlay and statistical significance
ggplot(metadata, aes(x = group_label, y = productive_templates, color = group_label)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(x = "", y = "Number of Productive TCRs") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("Short Interval", "No NACT")),
            map_signif_level = TRUE, color = "black")

# Create jitterplot with mean overlay and statistical significance
ggplot(metadata, aes(x = group_label, y = aa_richness, color = group_label)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(x = "", y = "Unique Amino Acid Sequences") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("No NACT", "Short Interval")), 
                                 #c("Short Interval", "Long Interval"),
                                 #c("No NACT", "Long Interval")),
              map_signif_level = TRUE, color = "black")

# Create jitterplot with mean overlay and statistical significance
ggplot(test_plot, aes(x = stringent, y = Distinct_IDs, color = stringent)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_boxplot(width = 0.2) +
  scale_color_manual(values = c("#9AC9E3", "#EFC3E6")) +
  labs(x = "", y = "Patient Repertoires Per Cluster") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Broadly Expressed", "Cancer Specific")) +
  geom_signif(comparisons = list(c("Broadly Expressed", "Cancer Specific")), 
              #c("Short Interval", "Long Interval"),
              #c("No NACT", "Long Interval")),
              map_signif_level = TRUE, color = "black") +
  scale_y_log10()

# Logo plots
ggseqlogo(cluster6063, method = 'prob' )
ggseqlogo(cluster3454, method = 'prob' )
ggseqlogo(cluster5380, method = 'prob' )
ggseqlogo(allSpecificClusters[allSpecificClusters %>% nchar() == 12], method = 'prob' )

