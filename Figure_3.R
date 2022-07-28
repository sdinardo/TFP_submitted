#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.

#  2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.  

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.



# Load packages and set working directory
library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(viridis)

setwd("~/Desktop/Lab/FCA/Follow-Up/Code")

#### CLUSTER WORK ####
# Import (almost complete) Seurat object created by Gabby Vida and Steve DiNardo
# This will be edited using the below code to create the final object you have
#testisv1 <- readRDS(file = "FCAloomToSeurat2TFP_Annotations.rds")

# Editing counts data so that the raw counts are now log-transformed (log2(x + 1)) for plotting
# This set of steps required more memory than my computer had, and was performed on the cluster
#counts.data <- GetAssayData(testisv1[["RNA"]], slot = "counts")
#log.counts <- Matrix(as.matrix(x = log2(counts.data + 1)), sparse = TRUE)
# Store this log-transformed counts data into a new "log.counts" assay
#log.counts_assay <- CreateAssayObject(counts = log.counts)
# Add the "log.counts" assay to the existing object
# This simultaneously creates a Key (logcounts_) that can be used to refer to the
# log.counts data when plotting, so as to override the default assay
#testisv1[["log.counts"]] <- log.counts_assay

# Now exporting the new testisv1 Seurat object so our final file includes the log-transformed counts
#saveRDS(testisv1, file = "FCAloomToSeurat2TFP_Annotations_Update.rds")

# This updated object is the final object available in "Input Files". 
# It is named "FCAloomToSeurat2TFP_Annotations.rds. You will not need
# to repeat the above code in the "CLUSTER WORK" section, and you will
# be able to find log-transformed counts in the "log.counts" assay.
# As clarified below, plotting can be done by switching default assays
# or using the Key "logcounts_" prior to the feature symbol.
#### END CLUSTER WORK ####

# Import Seurat object created by Gabby Vida and Steve DiNardo
testisv1 <- readRDS(file = "FCAloomToSeurat2TFP_Annotations.rds")

#### USING SEURAT OBJECT ####

# FIGURE 3A
genes_plot <- FeaturePlot(testisv1, features = "nFeature_RNA", pt.size = 1.0, cols = c("lightgrey", "navyblue"), order = TRUE) + ggtitle("genes per cell") + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15))
ggsave("genes_per_cell.png", plot = genes_plot, width = 20, height = 20, units = "cm")

# FIGURE 3B
umi_plot <- FeaturePlot(testisv1, features = "nCount_RNA", pt.size = 1.0, cols = c("lightgrey", "navyblue"), order = TRUE) + ggtitle("UMIs per cell") + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15))
ggsave("umis_per_cell.png", plot = umi_plot, width = 20, height = 20, units = "cm")

# FIGURE 3C
# This code was written by Gabby Vida
# Again, some of the below code was performed in the cluster because subsetting requires a lot of memory
germline <- subset(testisv1, idents = c("25", "22", "5", "78", "40", "97", "41", "51", "73", "35", "33", "48", "105", "64", "45", "13", "21", "56", "100", "3", "106", "16", "55", "63", "39", "27", "83", "15", "0", "28", "26", "24", "18", "50", "17", "96", "12", "30", "37", "75", "66", "10", "46"))

# Average number of UMIs and genes per clusters
# Extract 3 columns of meta.data
genes_UMI_byCluster <- germline@meta.data %>%
  group_by(leiden_res_6.0, nCount_RNA, nFeature_RNA)

# Calculate means
testPlot <- genes_UMI_byCluster %>%
  group_by(leiden_res_6.0) %>%
  summarise(aveUMI = mean(nCount_RNA), aveGenes = mean(nFeature_RNA))

# Make the cluster column a character vector
testPlot$leiden_res_6.0 <- as.character(testPlot$leiden_res_6.0)
# Then turn it back into a factor with the levels in the correct order
testPlot$leiden_res_6.0 <- factor(testPlot$leiden_res_6.0, levels = c("25", "22", "5", "78", "40", "97", "41", "51", "73", "35", "33", "48", "105", "64", "45", "13", "21", "56", "100", "3", "106", "16", "55", "63", "39", "27", "83", "15", "0", "28", "26", "24", "18", "50", "17", "96", "12", "30", "37", "75", "66", "10", "46"))

# Plot with labels
germline_umi_plot <- ggplot(testPlot, aes(leiden_res_6.0, aveUMI)) +
  geom_point(aes(size = aveGenes)) +
  labs(title = "Average UMI per cluster during germline progression", subtitle = "(see key for aveGenes per cluster)") +
  ylab("Average number of UMIs per Cluster") +
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.93, 0.6), axis.text = element_text(size = 18)) +
  scale_size(range = c(1, 9))

ggsave("germline_umi_plot.png", plot = germline_umi_plot, width = 40, height = 15, units = "cm")

# FIGURE 3D-E
# These plots were taken directly from ASAP (https://asap.epfl.ch/projects/ASAP24)
# This is also available at https://flycellatlas.org/ -> Testis -> ASAP
# To generate, select "Visualization" -> "Embedding_HVG_UMAP" -> "Controls: Coloring: Continuous"
# Then, "Data Type: Custom gene set" -> "Categorical gene metadata: _Chromosomes" -> "Category: X or Y"
# To use the graph, save SVG directly from the browser

# FIGURE 3F
# Note that this code uses an alternate assay of the Seurat object called "log.counts", one containing log-transformed raw counts
# This plot can be created by switching default assay: DefaultAssay(testisv1) <- "log.counts"
# or by using the Key (used below) to override the default "RNA" assay in specifying the feature argument
Mhc_plot <- FeaturePlot(testisv1, features = "logcounts_Mhc", max.cutoff = "q99", pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15)) + scale_color_viridis(option = "plasma")
ggsave("Mhc.png", plot = Mhc_plot, width = 20, height = 20, units = "cm")

# FIGURE 3G
# Note that this code uses an alternate assay of the Seurat object called "log.counts", one containing log-transformed raw counts
# This plot can be created by switching default assay: DefaultAssay(testisv1) <- "log.counts"
# or by using the Key (used below) to override the default "RNA" assay in specifying the feature argument
Hml_plot <- FeaturePlot(testisv1, features = "logcounts_Hml", max.cutoff = "q99", pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15)) + scale_color_viridis(option = "plasma")
ggsave("Hml.png", plot = Hml_plot, width = 20, height = 20, units = "cm")

# FIGURE 3I
# Use Seurat Object "Idents" of the (default) "leiden", the Leiden 0.4 dataset
Idents(testisv1) <- "leiden"

# Subset the object so that only the relevant cell types (clusters) are taken:
# a subset of the mature spermatocytes (2), muscle cells (4), hemocyte (13),
# neuron (28), trachea (20), mid-late cyst cells (8), hub cells (27)
testis_subset <- subset(testisv1, idents = c("2", "4", "13", "28", "20", "8", "27"))
levels(testis_subset) <- c("2", "4", "13", "28", "20", "8", "27")

# Plot features of interest such that pairs of markers for each cell type are together
somatic_expression <- DotPlot(testis_subset, features = c("CycB", "fzo", "tup", "upd1", "tbc", "eya", "trh", "grh", "sif", "Syt1", "pnr", "Hml", "CG33557", "Mhc"), cols = c("lightgrey", "navyblue"), dot.scale = 10) +
  coord_flip() +
  theme(text = element_text(size = 26), axis.text = element_text(size = 26), axis.title.x = element_text(margin = margin(t = 80)), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15))
# x-axis was hand-labeled in Adobe Illustrator

ggsave("somatic_expression_plot.png", plot = somatic_expression, width = 28, height = 25, units = "cm")

# FIGURE 3J
# Note that this code uses an alternate assay of the Seurat object called "log.counts", one containing log-transformed raw counts

# List of TFs is available in "Input Files" -> "TFs_list_500.txt"
# For each of the 500 TFs, first decide: is it exclusively in mature spermatocytes, both in other cell types and mature spermatocytes, or just non-germline tissues (based on the Li et al., 2022: Table S3)
# Ended up with 18 genes exclusively in mature spermatocytes, 127 in both spermatocytes and other cell types (germline or otherwise), and 355 not expressed in the germline
# 4 of the 355 genes in other cell types were not found in the Testis snRNA-seq data: sna, nerfin-1, CG7786, gem and were left out of the analysis (see below)

# First, we set a threshold of 0.2 (difference in average log2(counts + 1) for a gene in cluster 2 and cluster 3) as a boundary for upregulation relative to spermatogonia and early spermatocytes
# This threshold was determined to be about half of our set of empirical examples, which tend to be more obvious - a minimum of ~0.5 (Mhc, Hml, Syt1, grh, eya, upd1)
# Use "counts" slot to avoid exponentiation before averaging, which in this case tended to highlight noise
Idents(testisv1) <- "leiden"
avg_per_gene_examples <- data.frame(AverageExpression(testisv1, assays = "log.counts", slot = "counts", features = c("eya", "upd1", "Syt1", "grh", "Mhc", "Hml")))

# Compute the average expression (log2(counts + 1)) of each TF in each leiden cluster, then export the table
# Subtract average in spermatogonia and early spermatocytes (cluster 3, yellow in table) from average in mature spermatocytes (cluster 2, orange in table)
# See table for the cutoff at 0.2 - qualifying genes in green, non-qualifying genes in blue
TF_table <- read.table("TFs_list_500.txt", header = TRUE)
TF_list <- c(TF_table$GENE)
# 4 of the 500 TFs on this list were not found in this dataset, so I continued with the 496 other TFs
# The four TFs are: sna, nerfin-1, CG7786, gem
# Use "counts" slot to avoid exponentiation before averaging, which in this case tended to highlight noise
avg_per_gene <- data.frame(AverageExpression(testisv1, assays = "log.counts", slot = "counts", features = TF_list))
write.table(avg_per_gene, file = "TFs_avg_expression_leiden0.4.tsv", sep = "\t")
# This is Figure 3 - source data 1, which was moved to an .xlsx and sorted by a column created in Excel (Cluster 2 - Cluster 3)
