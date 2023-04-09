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
load("../data/processed/aa_table.RData")

# Search databases
db_search <- LymphoSeq2::searchDB(study_table = aa_table, dbname = "all", chain = "trb")
matched <- db_search %>% filter(antigen != "NA")
matched <- matched %>% merge(metadata %>% select(c("sample_name", "group_label")), 
                             by.y = "sample_name", by.x = "repertoire_id")

test <- matched %>% group_by(group_label, repertoire_id, pathology) %>% 
  summarize(count = sum(duplicate_count)) %>% arrange(-count)
test$pathology <- replace(test$pathology, 
                          test$pathology == "Human herpesvirus 4 (Epstein Barr virus)", 
                          "Epstein Barr virus (EBV)")
test$pathology <- gsub(pattern = fixed("Inflammatory bowel disease \\(IBD\\)\xa0"), 
                       replacement = "IBD", x = test$pathology)


ggplot(data = test, aes(x = group_label, y = count, fill = pathology)) +
  geom_bar(stat = "identity") +
  xlab("Group") +
  ylab("Antigen Specificity") +
  labs(x = "", y = "# Antigen Specific TCRs", fill = "Pathology") +
  theme(legend.text = element_text(size = 5), legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.5, "cm")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))


cancers <- c("Lung Cancer", "Melanoma", "Merkel cell carcinoma", "Neoantigen", 
             "Tumor associated antigen (TAA)", "Hepatocellular Carcinoma")

test <- matched %>% filter(pathology %in% cancers) %>% 
  group_by(group_label, repertoire_id, pathology) %>% 
  summarize(count = sum(duplicate_count)) %>% arrange(-count)

ggplot(data = test, aes(x = group_label, y = count, fill = pathology)) +
  geom_bar(stat = "identity") +
  xlab("Group") +
  ylab("Antigen Specificity") +
  labs(x = "", y = "# Antigen Specific TCRs", fill = "Pathology") +
  theme(legend.text = element_text(size = 5), legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.5, "cm")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))


test <- matched %>% filter(pathology %in% cancers) %>% group_by(repertoire_id, group_label) %>%
  summarize(n=sum(duplicate_count)) 
test$group_label <- factor(test$group_label,      # Reordering group factor levels
                         levels = c("No NACT", "Short Interval", "Long Interval"))

ggplot(data = test, aes(x = repertoire_id, y = n)) +
  geom_bar(stat = "identity") +
  facet_grid(~group_label, scales = "free_x", space = "free_x") +
  xlab("Group") +
  ylab("Antigen Specificity") +
  labs(x = "Samples", y = "# Antigen Specific TCRs", fill = "Pathology") +
  theme(legend.text = element_text(size = 5), legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.5, "cm")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
            