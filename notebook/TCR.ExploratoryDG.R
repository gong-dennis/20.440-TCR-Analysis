library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

##### Set up #####

# Read in TCR sequences
study_table <- LymphoSeq2::readImmunoSeq("../data/raw/ChemoProjTCRs/") %>% 
  topSeqs(top = 1000)
save(study_table, file = "../data/processed/study_table.RData")

# Read in metadata
metadata <- read.csv("../data/processed/precalculated.tsv", sep = "\t")

# Calculate summary statistics
summary_table <- LymphoSeq2::clonality(study_table)
save(summary_table, file = "../data/processed/summary_table.RData")

sample_names <- study_table %>% pull(repertoire_id) %>% unique()

# Extract productive sequences
aa_table <- LymphoSeq2::productiveSeq(study_table = study_table, 
                                      aggregate = "junction_aa", prevalence = TRUE)
save(aa_table, file = "../data/processed/aa_table.RData")

nuc_table <- LymphoSeq2::productiveSeq(study_table = study_table, aggregate = "junction", 
                                       prevalence = FALSE)
save(nuc_table, file = "../data/processed/nuc_table.RData")


samples <- aa_table %>% 
  dplyr::pull(repertoire_id) %>% unique()
LymphoSeq2::lorenzCurve(repertoire_ids = samples, study_table = aa_table)

LymphoSeq2::plotRarefactionCurve(study_table = aa_table)

##### Table 1 and Kaplan Meier #####

# Create a list of categorical variables
categorical_variables <- c("sex", "biopsy_location")

# Create a list of continuous variables
continuous_variables <- c("age", "total_templates", "survival_months", "nact_intvl")

# Create a stratified variable (e.g., treatment group)
strata <- "group_label"

table1 <- CreateTableOne(vars = continuous_variables,
                         strata = strata,
                         data = metadata)

# Print the table
print(table1, showAllLevels = TRUE, nonnormal = continuous_variables, exact = "stage", quote = FALSE, noSpaces = TRUE)

patients_data <- metadata %>% mutate(status = 1)
surv_object <- with(patients_data, Surv(survival_months, status))

# Fit the Kaplan-Meier model
km_fit <- survfit(surv_object ~ group_label, data = patients_data)

# Define custom theme for the risk table
custom_risk_table_theme <- theme(panel.background = element_blank(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())

# Define custom theme for the main plot
custom_main_plot_theme <- theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Plot the Kaplan-Meier curve using ggplot2
km_plot <- ggsurvplot(km_fit,
                      data = patients_data,
                      risk.table = FALSE,
                      pval = FALSE,
                      conf.int = FALSE,
                      legend.labs = c("Long NACT", "No NACT", "Short NACT"),
                      xlab = "Months After Surgery",
                      ylab = "Survival probability",
                      ggtheme = custom_main_plot_theme)
print(km_plot)


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


# Recurrently present TCRs

top_freq <- LymphoSeq2::topFreq(productive_table = aa_table, frequency = 0.001)
unique_seqs <- LymphoSeq2::uniqueSeqs(productive_table = aa_table)
sequence_matrix <- LymphoSeq2::seqMatrix(amino_table = aa_table, sequences = unique_seqs$junction_aa)
top_freq_matrix <- dplyr::full_join(top_freq, sequence_matrix)
top_freq_matrix %>% arrange(-numberSamples)

write.csv(top_freq_matrix %>% filter(numberSamples > 2), "../data/analysis/recurrentTCRs.csv")


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

NoNACT.IDs <- study_table %>% filter(group_label == "No NACT") %>% 
  pull(repertoire_id) %>% unique()
ShortNACT.IDs <- study_table %>% filter(group_label == "Short Interval") %>% 
  pull(repertoire_id) %>% unique()
LongNACT.IDs <- study_table %>% filter(group_label == "Long Interval") %>% 
  pull(repertoire_id) %>% unique()

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


ggplot(df_long, aes(x=x, y=value, color=variable)) + 
  geom_line() +
  labs(x="Sum of TCR Clusters", y="Cumulative Distribution Function", color="Variable") +
  theme_minimal() +
  theme(panel.grid=element_blank(), legend.position="top", 
        legend.text=element_text(size=7), legend.margin=margin(t=4)) +
  guides(color=guide_legend(title=NULL)) +# Remove legend title
  scale_color_manual(labels=c("No NACT", "Short Interval", "Long Interval"), values = c("#E69F00", "#56B4E9", "#009E73"))
  
ks.test(avgCDF_Short, avgCDF_Long, alternative = "greater")
ks.test(avgCDF_Long, avgCDF_No, alternative = "greater")

lapply(ShortNACT.IDs, getCancerClusters)

plot(avgCDF_Long, col = "blue", main = "Empirical CDFs")
lines(avgCDF_No, col = "red")

