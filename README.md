# Overview

This repo contains reproducible code used for analysis of TCRB sequencing data from patients with colorectal liver metastasis (CLM) treated with long or short doses of chemotherapy compared to non-treated patients [1]. The data are deposited on the immuneACCESS database at the following [link](https://clients.adaptivebiotech.com/pub/ad7a2d37-a0bc-4d88-813e-6dd7d762a65b) with DOI: 10.21417/EH2023GS. The goal of this project is to understand the role of chemotherapy in reshaping the adaptive immune response in cancer.

This repo was created by Dennis Gong and Thomas Usherwood as part of the Biological Networks class at MIT (20.440). Please direct questions to dgong3@mit.edu and thomasu@mit.edu

# Data

The data was generated using the [immunoSEQ hsTCRB Kit](https://www.immunoseq.com/) with samples obtained by Høye et al., with an abstract published in Cancer Research, accessible [here](https://aacrjournals.org/cancerres/article/82/12_Supplement/1346/699749/Abstract-1346-T-cell-receptor-repertoire). The dataset contains 92 samples, with 35 patients receiving neoadjuvant chemotherapy (NACT) for short interval and 15 patients receiving long interval NACT. An additional 35 patients did not receive NACT. All repertoires are stored in Adaptive ImmunoSEQ format. 

# Folder Structure

```
20.440-TCR-Analysis/
|__ README.md					<- this file
|__ driver.sh 					<- runs all
|__ src/ 						<- contains all scripts for generating results
	|__ data/ 					<- scripts for cleaning data
	|__ analysis/ 				<- scripts for producing results
	|__ visualization/ 			<- scripts for plotting results
	|__ util/ 					<- commonly reused scripts
|__ data/						<- contains all data used for project
	|__ raw/					<- contains raw data
		|__ ChemoProjTCRs/		<- contains all ImmunoSEQ formatted repertoires
	|__ processed/				<- contains any processed datasets
	|__ analysis/				<- contains data ready for visualization
		|__ SampleOverview.tsv	<- overview file with population summary statistics
|__ notebook/					<- exploratory notebooks
	|__ TCR.ExploratoryDG.r		<- Dennis's exploratory notebook
|__ fig/ 						<- contains figures
	|__ main_fig/				<- contains figures for main
	|__ tables/					<- contains all tables
	|__ supp_fig/				<- contains figures for supplement
```

# Installation

The driver file for the analysis in this pipeline is titled "driver.sh". The entire pipeline can be ran using: 
`./driver.sh`

Package dependencies can be found in "requirements.txt"

# References

1. Høye, E. et al. Abstract 1346: T cell receptor repertoire sequencing reveals chemotherapy-driven clonal expansion in colorectal liver metastases. Cancer Research 82, 1346 (2022).

