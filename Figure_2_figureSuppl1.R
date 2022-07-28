#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.

#  2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.  

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.


# Load packages and set working directory
library(Seurat)
library(ggplot2)

setwd("~/Desktop/Lab/FCA/Follow-Up/Code")

# Import Seurat object created by Gabby Vida and Steve DiNardo
testisv1 <- readRDS(file = "FCAloomToSeurat2TFP_Annotations.rds")

#### USING SEURAT OBJECT ####

# This code was written by Gabby Vida and Steve DiNardo

# FIGURE 2 S1A
full_6.0 <- DimPlot(testisv1, label = FALSE, pt.size = 2.0)
# Labeling was performed by hand in Adobe Illustrator
ggsave("full6.0_clusters.png", plot = full_6.0, width = 65, height = 60, units = "cm")

# For later use to produce the main figure counterpart to S2C (FIGURE 2L)
# Subset the Leiden 6.0 dataset to create an object with only the germline clusters
germline6.0 <- subset(testisv1, idents = c("25", "22", "5", "78", "40", "97", "41", "51", "73", "35", "33", "48", "105", "64", "45", "13", "21", "56", "100", "3", "106", "16", "55", "63", "39", "27", "83", "15", "0", "28", "26", "24", "18", "50", "17", "96", "12", "30", "37", "75", "66", "10", "46"))
# To set the order of clusters:
levels(germline6.0) <- c("25", "22", "5", "78", "40", "97", "41", "51", "73", "35", "33", "48", "105", "64", "45", "13", "21", "56", "100", "3", "106", "16", "55", "63", "39", "27", "83", "15", "0", "28", "26", "24", "18", "50", "17", "96", "12", "30", "37", "75", "66", "10", "46")

# FIGURE 2 S1B
# Reset Seurat Object "Idents" to be that of the Leiden 0.8 dataset
Idents(testisv1) <- "leiden_res_0.8"
full_0.8 <- DimPlot(testisv1, label = TRUE, pt.size = 2.0, label.size = 24)
ggsave("full0.8_clusters.png", plot = full_0.8, width = 58, height = 60, units = "cm")

# For later use to produce FIGURE S1C
# Subset the Leiden 0.8 dataset to create an object with only the germline clusters
germline0.8 <- subset(testisv1, idents = c("18", "5", "8", "9", "11", "6", "7", "0", "3"))
# To set the order of clusters:
levels(germline0.8) <- c("18", "5", "8", "9", "11", "6", "7", "0", "3")

# FIGURE 2 S1C
# Make a vector of features (genes) of interest
germline_markers <- c("p-cup", "Dic61B", "fzo", "CycB", "sa", "can", "aly", "Rbp4", "kmg", "aub", "bam", "zpg", "stg")

# Plot germline markers at lower (Leiden 0.8) resolution
germline0.8_markers <- DotPlot(germline0.8, features = germline_markers, cols = c("lightgrey", "navyblue"), dot.scale = 6) + coord_flip() +
  labs(title = "Expression by Leiden 0.8 cluster as germline progresses") +
  xlab("Germline Markers") +
  ylab("") + 
  theme(plot.title = element_text(size=15), axis.text.x = element_text(size = 13, angle = 45, hjust = 1), axis.title.y = element_text(size = 13))
ggsave("germline_markers_0.8.png", plot = germline0.8_markers, width = 24, height = 11, units = "cm")

# For FIGURE 2L, the same code to plot germline markers at Leiden 6.0 resolution
germline6.0_markers <- DotPlot(germline6.0, features = germline_markers, cols = c("lightgrey", "navyblue"), dot.scale = 6) + coord_flip() +
  labs(title = "Expression by cluster as germline progresses") +
  xlab("Germline Markers") +
  ylab("") + 
  theme(plot.title = element_text(size=15), axis.text.x = element_text(size = 13, angle = 45, hjust = 1), axis.title.y = element_text(size = 13))
ggsave("germline_markers_6.0.png", plot = germline6.0_markers, width = 34, height = 11, units = "cm")
