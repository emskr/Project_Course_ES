# Pancreatic cancer cell line GEM resistance analysis---------------------------

# Author: Emilis Skruodys
# Email: emisk975@student.liu.se

# Course code: 8BKG36
# Group: 7
# Dataset: 3

# Date: 2024-12-17

# Script name: Pancreatic cancer cell line GEM resistance analysis

################################################################################

# Working directory and libraries ##############################################
setwd("~/Desktop/Bioinformatics of Big Data/Working Directory")
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)

set.seed(1000)

# Script #######################################################################

# 1. Initial data preparation---------------------------------------------------

# Reading data
count_raw <- as.data.frame(read.table("edgeR_normcounts_count.tabular",
                                      header = TRUE))
count_annotations <- as.data.frame(read.table("edgeR_normcounts_annotations.tabular", 
                                              header = TRUE))

# Changing annotations column name for merging
colnames(count_annotations)[1] <- "GeneID"

# Merging data
count_merged <- merge.data.frame(count_raw, count_annotations, by = "GeneID")

# Reordering merged data frame columns
count_merged <- count_merged[, c(1, 8, 5, 6, 7, 2, 3, 4)]

# Renaming columns
colnames(count_merged) <- c("GeneID", "Gene", 
                            "Sensitive1", "Sensitive2", "Sensitive3", 
                            "Resistant1", "Resistant2", "Resistant3")

# NA value cleaning
count_na_clean <- na.omit(count_merged)

# Checking for duplicate genes and computing mean for duplicates 
# (Teachers recommendation and ChatGPT)
count_dupes_removed <- aggregate(. ~ Gene, data = count_na_clean, FUN = mean)

# Renaming rows
rownames(count_dupes_removed) <- count_dupes_removed[, 1]

# Creating final count data frame
count_final <- subset(count_dupes_removed, select = -c(1,2))

# Creating a Gene ID data frame
count_geneID <- subset(count_dupes_removed, select = -c(1, 3:8))

# 2. Metadata file generation---------------------------------------------------

metadata <- data.frame(
  SRA_accession_number = c("SRR10416711", "SRR10416712", "SRR10416713",
                           "SRR10416714", "SRR10416715", "SRR10416716"),
  Organism = rep("Homo sapiens", 6),
  Sample_cell_line = c("BxPC-31", "BxPC-32", "BxPC-33",
                       "BxPC-3-GR1", "BxPC-3-GR2", "BxPC-3-GR3"),
  Resistance = c(rep("Sensitive", 3), rep("Resistant", 3)),
  Sample_name = c("Sensitive1", "Sensitive2", "Sensitive3",
                 "Resistant1", "Resistant2", "Resistant3"),
  Sequencing = rep("RNA-Seq paired end", 6)
)

# 3. Unsupervised clustering----------------------------------------------------

# Applying PCA on the transposed data
output_pca <- prcomp(t(count_final))

# Calculating the variance of each PC
var_explained = output_pca$sdev^2 / sum(output_pca$sdev^2)
var_explained <- round(var_explained * 100, 2)

# Showing the first 6 PCs
var_explained_df <- data.frame(PC = paste0("PC", 1:6),
                               var_explained = var_explained[1:6])

# Extracting PCs 1 and 2
output_pca <- as.data.frame(output_pca$x[, 1:2])

# Adding metadata variables to the PCA output
output_pca$Resistance <- metadata$Resistance
output_pca$Sample_name <- metadata$Sample_name

# Creating a folder for all plots
dir.create("Plots")

# Plotting samples by their PC1 and PC2 (with ChatGPT)
ggplot(output_pca, aes(x = PC1, 
                       y = PC2, 
                       color = Resistance)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample_name), 
            vjust = -1, 
            hjust = 0.5, 
            size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("PCA plot") +
  xlab(paste0("PC1 (", var_explained[1], "%)")) +
  ylab(paste0("PC2 (", var_explained[2], "%)"))

# Saving plot as .png
ggsave("1_PCA.png",
  plot = last_plot(),
  path = "Plots/",
  width = 7,
  height = 5)

# 4. Statistical testing--------------------------------------------------------

# Creating a P value data frame
pvalue <- data.frame(p_value = rep(0, nrow(count_final)))
rownames(pvalue) <- rownames(count_final)

# Performing student's t-test
for (a in 1:nrow(count_final)){
  pvalue$p_value[a] <- t.test(count_final[a,1:3],
                              count_final[a,4:6],
                              alternative = "two.sided")$p.value
}

# Significant P value filtering
sig_pval <- pvalue[pvalue$p_value < 0.05,,drop = FALSE]

# Subsetting significant genes after t-test
names_sig_genes <- rownames(sig_pval)
sig_genes <- subset(count_final, rownames(count_final) %in% names_sig_genes)

# 5. Correction for multiple tests----------------------------------------------

# Performing BH correction
pvalue_bh <- data.frame(p_value = p.adjust(sig_pval$p_value, method = "BH")) # With ChatGPT
rownames(pvalue_bh) <- rownames(sig_genes)

# Changing column name
colnames(pvalue_bh) <- "p_value_adj"

# Filtering significant genes after BH correction
sig_pval_bh <- pvalue_bh[pvalue_bh$p_value < 0.05,,drop = FALSE]

# 6. Calculating log fold change------------------------------------------------

# Creating empty data frame for filtered count means
means_fc <- data.frame(sensitive_means = rep(0,nrow(sig_genes)),
                       resistant_means = rep(0,nrow(sig_genes)))

# Adding mean count values to the empty mean data frame 
for (b in 1:nrow(sig_genes)){
  means_fc$sensitive_means[b] <- rowMeans(sig_genes[b,1:3])
  means_fc$resistant_means[b] <- rowMeans(sig_genes[b,4:6])
}

# Calculating log fold change for each gene (count data is already log adjusted by Galaxy)
for (c in 1:nrow(means_fc)){
  means_fc$logFC[c] <-  means_fc[c,2] - means_fc[c,1]
}

# Adding gene names
rownames(means_fc) <- rownames(sig_genes)

# Creating a data frame with GeneID
sig_genes_ID <- sig_genes
sig_genes_ID$GeneID <- count_geneID[row.names(sig_genes_ID), "GeneID"]

# Filtering means_fc_filter by log fold change threshold
means_fc_filter <- means_fc[!(means_fc$logFC > -1 & means_fc$logFC < 1), ]

# Removing genes with an insignificant log fold change
sig_genes_final <- sig_genes[rownames(sig_genes) %in% rownames(means_fc_filter), ]

# 7. Creating unified data frame------------------------------------------------

# Creating unified data frame
data_final <- means_fc_filter

# Subsetting rows in pvalue_bh to match row names in data_final
pvalue_bh_filtered <- pvalue_bh[rownames(data_final),, drop = FALSE]

# Combining the data frames
data_final <- cbind(data_final, pvalue_bh_filtered)

# Subsetting rows in count_final to match row names in data_final
sig_genes_filtered <- sig_genes[rownames(data_final),, drop = FALSE]

# Combining the data frames by columns
data_final <- cbind(data_final, sig_genes_filtered)

# Subsetting rows in count_final to match row names in data_final
count_geneID_filtered <- count_geneID[rownames(data_final),, drop = FALSE]

# Combining the data frames by columns
data_final <- cbind(data_final, count_geneID_filtered)

# Reordering data
data_final <- data_final[, c("Sensitive1", "Sensitive2", "Sensitive3", 
                             "Resistant1", "Resistant2", "Resistant3",
                             "sensitive_means", "resistant_means",
                             "p_value_adj", "logFC", "GeneID")]

# 8. Visualization--------------------------------------------------------------

# 8.1. Volcano plot----

# Creating an empty data frame for volcano plot
volcano_plot_data <- data.frame(sensitive_means = rep(0,nrow(count_final)),
                                resistant_means = rep(0,nrow(count_final)))

# Filling the volcano plot with unfiltered data
for (d in 1:nrow(count_final)){
  volcano_plot_data$sensitive_means[d] <- rowMeans(count_final[d,1:3])
  volcano_plot_data$resistant_means[d] <- rowMeans(count_final[d,4:6])
  
}

for (e in 1:nrow(volcano_plot_data)){
  volcano_plot_data$logFC[e] <-  volcano_plot_data[e,1] - volcano_plot_data[e,2]
}

# Adding non-adjusted P values
volcano_plot_data <- cbind(volcano_plot_data, pvalue)

# Renaming columns according to data
colnames(volcano_plot_data) <- c("sensitive_means", "resistant_means", 
                                 "logFC_unfiltered", "p_value_unfiltered")

# Adjusting p values using BH method
volcano_plot_data$p_value_unfiltered <- p.adjust(volcano_plot_data$p_value_unfiltered, 
                                                 method = "BH")

# Adding differential expression
volcano_plot_data$diff_expressed <- "No differential expression"

volcano_plot_data$diff_expressed[volcano_plot_data$logFC_unfiltered > 1 &
                                   volcano_plot_data$p_value_unfiltered < 0.05] <- "Up-regulated"

volcano_plot_data$diff_expressed[volcano_plot_data$logFC_unfiltered < -1 & 
                                   volcano_plot_data$p_value_unfiltered < 0.05] <- "Down-regulated"

# Plotting the volcano plot
ggplot(data = volcano_plot_data, aes(x = logFC_unfiltered, 
                                     y = -log10(p_value_unfiltered),
                                     col = diff_expressed)) +
  scale_color_manual(values = c("brown", "grey", "red")) + 
  ggtitle("Volcano Plot of Differential Gene Expression") +
  labs(y = "-Log10 (P value adj.)", 
       x = "Log2 fold change",
       colour = "Differential expression") +
  geom_point(size = 0.8) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed")+
  theme_classic()

# Saving plot as .png
ggsave("2_Volcano.png",
  plot = last_plot(),
  path = "Plots/",
  width = 10,
  height = 7)

# 8.2. Lollipop----
# (With Chat GPT)

# Extracting top 10 sensitive means
top10_sensitive <- data_final[order(-data_final$sensitive_means, 
                                    na.last = TRUE), ][1:10, ]
top10_sensitive$type <- "Sensitive"

# Extracting top 10 resistant means
top10_resistant <- data_final[order(-data_final$resistant_means, 
                                    na.last = TRUE), ][1:10, ]
top10_resistant$type <- "Resistant"

# Combining the two datasets
combined_data <- rbind(
  data.frame(Sample = row.names(top10_sensitive), 
             Value = top10_sensitive$sensitive_means, 
             Type = top10_sensitive$type),
  data.frame(Sample = row.names(top10_resistant), 
             Value = top10_resistant$resistant_means, 
             Type = top10_resistant$type)
)

# Creating the top lollipop plot
ggplot(combined_data, aes(x = reorder(Sample, Value), 
                          y = Value, 
                          color = Type)) +
  geom_segment(aes(xend = Sample, 
                   yend = 0),
               linewidth = 1) +
  geom_point(size = 4) +
  coord_flip() +
  labs(
    title = "Top 10 Sensitive and Resistant Gene Expression Means",
    x = "Genes",
    y = "Gene expression",
    color = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Saving plot as .png
ggsave("3_Top10_lollipop.png",
       plot = last_plot(),
       path = "Plots/",
       width = 8,
       height = 8)

# Extracting bottom 10 rows based on sensitive_means
bottom10_sensitive <- data_final[order(data_final$sensitive_means, 
                                       na.last = TRUE), ][1:10, ]
bottom10_sensitive$type <- "Sensitive"

# Extracting bottom 10 rows based on resistant_means
bottom10_resistant <- data_final[order(data_final$resistant_means, 
                                       na.last = TRUE), ][1:10, ]
bottom10_resistant$type <- "Resistant"

# Combining the two datasets
combined_bottom_data <- rbind(
  data.frame(Sample = row.names(bottom10_sensitive), 
             Value = bottom10_sensitive$sensitive_means, 
             Type = bottom10_sensitive$type),
  data.frame(Sample = row.names(bottom10_resistant), 
             Value = bottom10_resistant$resistant_means, 
             Type = bottom10_resistant$type)
)

# Creating the bottom lollipop plot
ggplot(combined_bottom_data, 
       aes(x = reorder(Sample, Value), 
           y = Value, 
           color = Type)) +
  geom_segment(aes(xend = Sample, yend = 0), linewidth = 1) +
  geom_point(size = 4) +
  coord_flip() +
  scale_y_continuous(
    breaks = seq(-12.5, 0, by = 2.5),  #
    limits = c(-12.5, 0)) +
  labs(title = "Bottom 10 Sensitive and Resistant Gene Expression Means",
       x = "Genes",
       y = "Gene expression",
       color = "Category") +
  theme_minimal() 

# Saving plot as .png
ggsave("4_Bottom10_lollipop.png",
       plot = last_plot(),
       path = "Plots/",
       width = 8,
       height = 8)

# 8.3. Heatmap----

# Creating ordered P value data frame
data_final_p_ordered <- as.data.frame(data_final[order(data_final$p_value_adj),,
                                                 drop = FALSE])
data_final_p_ordered <- subset(data_final_p_ordered, select =  -c(7:11))

# Subsetting to leave top 40 significant genes by P value
heatmap_data <- as.data.frame(data_final_p_ordered[1:40,,drop = FALSE])

# Saving plot as .png
png("Plots/5_Heatmap.png",
    width = 1300,
    height = 1600)

# Plotting the heatmap
pheatmap(heatmap_data,
         fontsize_row = 25,
         fontsize_col = 25,
         fontsize = 40,
         cellwidth = 130, 
         cellheight = 35,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         angle_col = 0,
         legend = FALSE,
         main = "Heatmap of Top 40 Significant Gene Expression"
)

dev.off()

# 8.4. Pathway enrichment analysis using KEGG----

# Creating file for enrichment analysis
enrichment_data <- subset(data_final, select = -c(1:9))

# Renaming rows as GeneID
row.names(enrichment_data) <- enrichment_data$GeneID

# Subsetting by differential expression
logFC_up <- subset(enrichment_data, enrichment_data$logFC > 1)

# Performing KEGG enrichment for up-regulated genes
kegg_up <- enrichKEGG(gene = logFC_up$GeneID, 
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.05, 
                        pAdjustMethod = "BH")

# Plotting KEGG enrichment analysis for up-regulated pathways
dotplot(kegg_up,
        showCategory = 10,
        title = "Top KEGG up-regulated enriched pathways")

# Saving plot as .png
ggsave("6_KEGG_up-regulated.png",
       plot = last_plot(),
       path = "Plots/",
       width = 8,
       height = 8)

# Performing KEGG enrichment for all genes
kegg_all <- enrichKEGG(gene = data_final$GeneID,
                       organism = "hsa",
                       keyType = "kegg",
                       pvalueCutoff = 0.05, 
                       pAdjustMethod = "BH")

# Plotting KEGG enrichment analysis for all gene pathways
dotplot(kegg_all,
        showCategory = 10,
        title = "Top KEGG enriched pathways")

# Saving plot as .png
ggsave("7_KEGG_all.png",
       plot = last_plot(),
       path = "Plots/",
       width = 8,
       height = 8)
