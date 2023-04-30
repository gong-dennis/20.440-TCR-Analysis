# Driver script used to perform analysis and generate figures

# 1. Preprocessing

Rscript src/analysis/SaveSummaryTables.R # Saves several essential tables for downstream analysis
Rscript src/analysis/CancerTCRSummaryStats.R # Creates "metadata" dataframe used for jitterplots

# 2. Generate Main figures & tables
# 	 Figures stored in fig/main_fig/
# 	 Tables stored in fig/tables/

Rscript src/analysis/CDF.R # Creates cumulative distribution function plot
Rscript src/visualization/StackedBarPrevalence.R # Generate cancer associated TCR stacked bar chart
Rscript src/analysis/ClusterSizes.R # Generate histogram of cluster sizes
Rscript src/visualization/LogoPlots.R # Generate histogram of cluster sizes

Rscript src/visualization/Fig1A.R # Generate histogram of cluster sizes
Rscript src/visualization/Fig1B.R # Generate histogram of cluster sizes
Rscript src/visualization/Fig1C.R # Generate histogram of cluster sizes

Rscript src/visualization/Fig2A.R # Generate histogram of cluster sizes
Rscript src/visualization/Fig2B.R # Generate histogram of cluster sizes
Rscript src/visualization/Fig2C.R # Generate histogram of cluster sizes
Rscript src/visualization/Fig4D.R # Generate histogram of cluster sizes

Rscript src/analysis/NetworkAnalysis.R # Generate histogram of cluster sizes

# 3. Generate Supplementary figures & tables
# 	 Stored in fig/supp_fig/

Rscript src/visualization/aaLength.R # Generate supplemental aa length plot
Rscript src/analysis/VGenePlots.R # Show V gene diversity across samples

python3 main.py
