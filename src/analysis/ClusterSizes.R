library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(forcats)

##### Clustering all samples together #####

allSampleCluster <- read.csv("data/processed/AllSampleCluster.csv")

stringentPlot <- allSampleCluster %>% group_by(cluster, stringentCancer) %>% summarize(n=n())
cancerPlot <- allSampleCluster %>% group_by(cluster, cancer) %>% summarize(n=n())

ggplot(cancerPlot, aes(x = n, fill = cancer)) +
  geom_histogram(binwidth = 5, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("#9AC9E3", "#EFC3E6"), name = "", 
                    labels = c("Broadly Expressed", "Cancer Associated")) +
  labs(x = "Number of TCRs per Cluster", y = "Frequency") +
  theme_classic() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(1, 10, 100, 1000),
                     labels = c(1, 10, 100, 1000)) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))
