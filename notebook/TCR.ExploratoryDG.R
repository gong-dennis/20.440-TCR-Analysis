library(tidyverse)
library(LymphoSeq2)
library(ggplot2)

setwd("~/Documents/GitHub/20.440-TCR-Analysis/notebook")

# Read in TCR sequences
study_table <- LymphoSeq2::readImmunoSeq("../data/raw/ChemoProjTCRs/") %>% 
  topSeqs(top = 1000)


summary_table <- LymphoSeq2::clonality(study_table)
summary_table

sample_names <- study_table %>% pull(repertoire_id) %>% unique()

# Extract productive sequences
aa_table <- LymphoSeq2::productiveSeq(study_table = study_table, 
                                      aggregate = "junction_aa", prevalence = TRUE)

nuc_table <- LymphoSeq2::productiveSeq(study_table = study_table, aggregate = "junction", 
                                       prevalence = FALSE)

# Clonal relatedness

LymphoSeq2::clonalRelatedness(study_table = study_table, edit_distance = 10)

LymphoSeq2::topSeqsPlot(study_table = aa_table, top = 10)

top_freq <- LymphoSeq2::topFreq(productive_table = aa_table, frequency = 0.001)
unique_seqs <- LymphoSeq2::uniqueSeqs(productive_table = aa_table)
sequence_matrix <- LymphoSeq2::seqMatrix(amino_table = aa_table, sequences = unique_seqs$junction_aa)
top_freq_matrix <- dplyr::full_join(top_freq, sequence_matrix)
top_freq_matrix


# Chord diagram VDJ
top_seqs <- LymphoSeq2::topSeqs(nuc_table, top = 1)
LymphoSeq2::chordDiagramVDJ(study_table = top_seqs, 
                            association = "VJ", 
                            colors = c("darkred", "navyblue"))


vGenes <- LymphoSeq2::geneFreq(nuc_table, locus = "V", family = TRUE)
multicolors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Set1")))(28)
ggplot2::ggplot(vGenes, aes(x = repertoire_id, y = gene_frequency, fill = gene_name)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::theme_minimal() + 
  ggplot2::scale_y_continuous(expand = c(0, 0)) + 
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::scale_fill_manual(values = multicolors) + 
  ggplot2::labs(y = "Frequency (%)", x = "", fill = "") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

# Search databases

db_search <- LymphoSeq2::searchDB(study_table = aa_table, dbname = "all", chain = "trb")
colnames(db_search)
