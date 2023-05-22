library(tidyverse)
library(LymphoSeq2)
library(ggplot2)
library(gridExtra)
library(tableone)
library(ggstance)
library(tidytext)
library(survminer)

# Read in metadata
metadata <- read.csv("data/processed/precalculated.tsv", sep = "\t")

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


# Generate Kaplan Meier Curve
patients_data <- metadata %>% mutate(status = 1)
surv_object <- with(patients_data, survival::Surv(survival_months, status))

# Fit the Kaplan-Meier model
km_fit <- survival::survfit(surv_object ~ group_label, data = patients_data)

# Define custom theme for the main plot
custom_main_plot_theme <- theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

# Plot the Kaplan-Meier curve using ggplot2
ggsurvplot(km_fit,
           data = patients_data,
           palette = c("#E69F00", "#56B4E9", "#009E73"),
           risk.table = FALSE,
           pval = TRUE,
           conf.int = FALSE,
           legend.labs = c("Long NACT", "No NACT", "Short NACT"),
           xlab = "Months After Surgery",
           ylab = "Survival probability",
           ggtheme = custom_main_plot_theme)

