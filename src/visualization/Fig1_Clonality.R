library(tidyverse)
library(ggplot2)

# From Adaptive summary file
study_statistics <- read.table("data/analysis/SampleOverview.tsv", 
                               header = TRUE, sep = '\t')

# Visualize histogram
p1 <- ggplot(study_statistics, aes(x = productive_simpson_clonality)) +
  geom_histogram()
print(p1)

# Save to pdf

savedir <- "fig/supp_fig/S1.1_prodSimpClon.pdf"

pdf(file = savedir,   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

print(p1)

dev.off()

