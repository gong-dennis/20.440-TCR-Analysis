library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(survminer)

load(file = "data/processed/study_table.RData")

# Read in metadata
metadata <- read.csv("data/processed/precalculated.tsv", sep = "\t")

##### Generate aa length heatmap #####

study_table <- merge(study_table, metadata %>% select(sample_name, group_label), by.x = "repertoire_id", by.y = "sample_name")
study_table <- study_table[is.finite(study_table$junction_aa_length), ]

# Bin the length measurements
data_binned <- study_table %>%
  mutate(bin = floor(junction_aa_length)) %>%
  group_by(repertoire_id, group_label, bin) %>%
  summarize(count = n(), .groups = "drop")
data_binned <- data_binned %>%
  group_by(repertoire_id) %>%
  mutate(total = sum(count), proportion = count / total) %>%
  ungroup()
data_binned$repertoire_id <- factor(data_binned$repertoire_id)
data_binned$group_label <- factor(data_binned$group_label, levels = c("No NACT", "Short Interval", "Long Interval"))
data_binned$bin <- factor(data_binned$bin)

# Create the aa length heatmap
ggplot(data_binned, aes(x = repertoire_id, y = bin, fill = proportion)) +
  geom_tile() +
  scale_fill_gradient(low = "grey", high = "steelblue", na.value = "white", n.breaks = 4) +
  facet_grid(~group_label, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Sample ID", y = "CDR3 AA Length", fill = "Proportion of sequences")

