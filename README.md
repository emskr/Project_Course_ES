# Project_Course_ES
This GitHub repository contains the R code and data used in the project of the Project Course: Bioinformatics of Big Data.

The data files in this repository are named 'edgeR_normcounts_count.tabular', which consist of the normalized and log-transformed gene read counts, and "edgeR_normcounts_annotations.tabular", which includes the gene IDs and their annotations (gene names) of the generated gene counts in the previous file.

These data files were generated using the web-based platform Galaxy (https://usegalaxy.eu/) by pre-processing the RNA-seq data generated in the study "Identification of chemoresistance-related mRNAs based on gemcitabine-resistant pancreatic cancer cell lines" by Zhou et al., 2019.

This project analyzes and compares transcriptomic data to identify transcriptional response changes leading to GEM resistance in the BxPC-3 pancreatic cell line.

The NIH SRA database SRR accession numbers used in this project were:
- SRR10416711 <br/>
- SRR10416712 <br/>
- SRR10416713 <br/>
- SRR10416714 <br/>
- SRR10416715 <br/>
- SRR10416716 <br/>

The pre-processing performed in Galaxy using the following tools:
- Accession of SRR files (using Faster Download and Extract Reads in FastQ) followed by quality control <br/>
- Quality control (using MultiQC) <br/>
- Adapter trimming (using fastp) followed by quality control <br/>
- Genome alignment (using HISAT2) followed by quality control <br/>
- Read quantification (using featureCounts) <br/>
- Differential expression analysis (using edgeR) <br/>

The code in R contains the following:
- Initial data preparation (read-in and data formatting) <br/>
- Generation of a metadata file with relevant information about the samples, including the SRA accession number, organism, sample cell line, GEM resistance status, sample names used in the script and generated plots, and sequencing method. <br/>
- Unsupervised clustering using principal component analysis <br/>
- Statistical testing using student's t-test to extract P values <br/>
- Calculation of log fold change <br/>
- Adjustment of P values using the Benjamini-Hochberg method <br/>
- Creation of a unified data frame containing generated data required in further steps <br/>
- Visualisation in different plots: <br/>
  - Volcano plot for visualising differential expression <br/>
  - Lollipop plots for comparing the most and least expressed genes <br/>
  - Heatmap for the comparison of the most signifcant gene expression <br/>
  - Pathway enrichment analysis plots using KEGG <br/>
