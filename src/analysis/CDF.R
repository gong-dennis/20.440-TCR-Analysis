library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggstance)
library(tidytext)

# Read in metadata
metadata <- read.csv("data/processed/precalculated.tsv", sep = "\t")

load(file = "data/processed/aa_table.RData")

##### Characterize cluster distributions in each group #####

getCancerClusters <- function(name) {
  clusters <- read.csv(paste0("data/processed/cancer_clusters/", name, ".csv")) 
  if ("cluster" %in% names(clusters)) {
    return(clusters %>% pull(cluster))
  } else {
    return(vector())
  }
}

getNumTCRsFromCluster <- function(name, clusters) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% 
    merge(TCRs, by = "junction_aa") %>% filter(cluster %in% clusters)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  return(table)
}

getCDF <- function(name) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% merge(TCRs)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  return(ecdf(table$n))
}

getTCRsPerCluster <- function(name) {
  TCRs <- read.csv(paste0("data/processed/clustcr_labels/", name, ".csv"))
  mergedAA <- as_tibble(aa_table) %>% filter(repertoire_id == name) %>% merge(TCRs)
  table <- mergedAA %>% group_by(cluster) %>% summarize(n = sum(duplicate_count)) %>%
    arrange(-n)
  return(table$n)
}

# Generate CDFs for each group

NoNACT.IDs <- metadata %>% filter(group_label == "No NACT") %>% 
  pull(sample_name) %>% unique()
ShortNACT.IDs <- metadata %>% filter(group_label == "Short Interval") %>% 
  pull(sample_name) %>% unique()
LongNACT.IDs <- metadata %>% filter(group_label == "Long Interval") %>% 
  pull(sample_name) %>% unique()


CDFs <- lapply(metadata$sample_name, getCDF)
names(CDFs) <- metadata$sample_name

x_values <- seq(0, 1000, length.out=100)
avgCDF_No <- sapply(x_values, function(x) mean(sapply(CDFs[NoNACT.IDs], function(f) f(x))))
avgCDF_Short <- sapply(x_values, function(x) mean(sapply(CDFs[ShortNACT.IDs], function(f) f(x))))
avgCDF_Long <- sapply(x_values, function(x) mean(sapply(CDFs[LongNACT.IDs], function(f) f(x))))

df_avg_cdf <- data.frame(x = x_values, `No NACT` = avgCDF_No, 
                         `Short Interval` = avgCDF_Short, `Long Interval` = avgCDF_Long)
df_long <- gather(df_avg_cdf, key="variable", value="value", -x)
df_long$variable <- factor(df_long$variable, levels=c("No.NACT", "Short.Interval", "Long.Interval"))

# Plot
ggplot(df_long, aes(x=x, y=value, color=variable)) + 
  geom_line() +
  labs(x="Cumulative TCR Clusters", y="Cumulative Probability", color="Variable") +
  theme_minimal() +
  theme(panel.grid=element_blank(), legend.position="top", 
        legend.text=element_text(size=7, face = "bold"), legend.margin=margin(t=4),
        axis.title = element_text(face = "bold")) +
  guides(color=guide_legend(title=NULL)) +# Remove legend title
  scale_color_manual(labels=c("No NACT", "Short Interval", "Long Interval"), values = c("#E69F00", "#56B4E9", "#009E73"))

# Test for statistical significance
ks.test(avgCDF_Short, avgCDF_Long, alternative = "greater")
ks.test(avgCDF_Long, avgCDF_No, alternative = "greater")

