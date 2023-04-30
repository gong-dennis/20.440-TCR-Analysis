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
ggplot(metadata, aes(x = survivor, y = propCancer, color = survivor)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  scale_color_manual(values = c("#C1E1C5", "#FF9F80")) +
  labs(x = "", y = "Proportion Cancer Associated TCRs") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face="bold")) +
  geom_signif(comparisons = list(c("Short Term", "Long Term")), 
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
