library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(forcats)

# Calculate prevalence of cancer associated TCRs and plot

caTCRPrevalenceTable <- function(name) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% merge(TCRs)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  TCRs <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  table$cancer <- table$cluster %in% TCRs$cluster
  table <- head(table, 20)
  table$prop <- table$n / sum(table$n)
  table$sample_name <- name
  return(arrange(table, -prop))
}

metadata <- read.csv("data/processed/precalculated.tsv", sep = "\t")
sample_names <- metadata$sample_name

hi <- lapply(sample_names, caTCRPrevalenceTable) %>% do.call(what = rbind)
hi <- hi %>% merge(select(metadata, c(sample_name, group_label)), by = "sample_name")
hi$cancer[hi$cancer == TRUE] <- "Cancer Associated"
hi$cancer[hi$cancer == FALSE] <- "Non-Cancer Associated"

my_colors <- c("#9AC9E3", "#EFC3E6")
ggplot(hi, aes(x = sample_name, y = prop,
               fill = fct_reorder(factor(cancer), prop),
               color = cancer)) +
  geom_col(position = "fill", colour="black", size = 0.2) +
  facet_grid(~ group_label, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "", y = "Proportion of 20 Top Clusters", fill = "") +
  theme_classic() + ylim(0,1) +
  theme(legend.position = "top") +
  scale_fill_manual(values = my_colors)
