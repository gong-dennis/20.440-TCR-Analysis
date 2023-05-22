library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(forcats)
library(ggsignif)

# Load in data
load(file = "data/processed/aa_table.RData")
allSampleCluster <- read.csv("data/processed/AllSampleCluster.csv")

# Match clusters with amino acids
test <- merge(allSampleCluster, aa_table, by = "junction_aa")

stringentCancerTCRs <- test %>% group_by(cluster) %>% summarize(stringentCancer = sum(stringentCancer == "Yes"),
                                                                stringentNormal = sum(stringentNormal == "Yes"),
                                                                count = sum(duplicate_count), 
                                                                repertoires = n_distinct(repertoire_id)) %>%
  filter(stringentCancer > stringentNormal) %>% arrange(-stringentCancer) %>%
  filter(stringentNormal == 0) %>% filter(repertoires > 1)

test_plot <- test %>% group_by(cluster, cancer) %>% summarize(Distinct_IDs = n_distinct(repertoire_id), n = n()) %>% arrange(-Distinct_IDs)

# Create a category column
test_plot$stringent <- "Broadly Expressed"
test_plot$stringent[test_plot$cluster %in% stringentCancerTCRs$cluster] <- "Cancer Specific"

# Create jitterplot with mean overlay and statistical significance
ggplot(test_plot, aes(x = stringent, y = Distinct_IDs, color = stringent)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_boxplot(width = 0.2) +
  scale_color_manual(values = c("#9AC9E3", "#EFC3E6")) +
  labs(x = "", y = "Patient Repertoires Per Cluster") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face="bold")) +
  scale_x_discrete(labels = c("Broadly Expressed", "Cancer Specific")) +
  geom_signif(comparisons = list(c("Broadly Expressed", "Cancer Specific")), 
              #c("Short Interval", "Long Interval"),
              #c("No NACT", "Long Interval")),
              map_signif_level = TRUE, color = "black") +
  scale_y_log10()

ggplot(test_plot, aes(x = Distinct_IDs, fill = stringent)) +
  geom_histogram(binwidth = 1, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("#9AC9E3", "#EFC3E6"), name = "", 
                    labels = c("Broadly Expressed", "Cancer Specific")) +
  labs(x = "Patient Repertoires Per Cluster", y = "Frequency") +
  theme_classic() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(1, 10, 100, 1000),
                     labels = c(1, 10, 100, 1000)) +

  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))

