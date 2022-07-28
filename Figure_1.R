# Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581). Code provided ‘as is’, under terms of the MIT license.

# 2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.

# Load packages and set working directory
library(Seurat)
library(ggplot2)

setwd("~/Desktop/Lab/FCA/Follow-Up/Code")

# Import Seurat object created by Gabby Vida and Steve DiNardo
testisv1 <- readRDS(file = "FCAloomToSeurat2TFP_Annotations.rds")

#### USING SEURAT OBJECT ####

# FIGURE 1B
cell_types_plot <- DimPlot(testisv1, cols = c("#8f00b3", "#3d1040", "#d580ff", "#fff700", "#fff700", 
                           "#fff700", "#fff700", "#fff700", "#fff700", "#0433ff", "#0433ff", 
                           "#fff700", "#149900", "#a1ff80", "#0433ff", "#cbb4ff", "#0433ff", 
                           "#0433ff", "#0433ff", "#accbe6", "#b23000", "#4b008c", "#0433ff", 
                           "#0433ff", "#0433ff", "#0433ff", "#0433ff", "#0433ff", "#0433ff", 
                           "#0433ff", "#0433ff", "#0433ff", "#fff700", "#fff700", "#0433ff", 
                           "#0433ff", "#4c4359", "#fff700", "#a6538a", "grey"), 
                           pt.size = 1.0, group.by = "TFP_Annotations", label = FALSE, order = FALSE) + NoLegend() + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"))
ggsave("cell_types.png", plot = cell_types_plot, width = 18, height = 20, units = "cm")

# FIGURE 1C
tj_vas_plot <- FeaturePlot(testisv1, features = c("tj", "vas"), cols = c("lightgrey", "#fde725ff", "#1502B9"), blend = TRUE, blend.threshold = 0, pt.size = 1.0, order = TRUE) & theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 20))
ggsave("tj_vas.png", plot = tj_vas_plot, width = 75, height = 20, units = "cm")

# FIGURE 1D
stg_plot <- FeaturePlot(testisv1, features = "stg", cols = c("lightgrey", "navyblue"), pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15))
ggsave("stg.png", plot = stg_plot, width = 20, height = 20, units = "cm")

# FIGURE 1E
esg_plot <- FeaturePlot(testisv1, features = "esg", cols = c("lightgrey", "navyblue"), pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 20))
ggsave("esg.png", plot = esg_plot, width = 20, height = 20, units = "cm")

# FIGURE 1F
annotations_plot <- DimPlot(testisv1, cols = c("#8f00b3", "#3d1040", "#d580ff", "#4c3213", "#005359", 
                            "#79eaf2", "#99737d", "#f23d55", "#33262a", "#7f0011", "#f26100", 
                            "#59a1b3", "#149900", "#a1ff80", "#583973", "#cbb4ff", "#79bf60", 
                            "#a65800", "#7ca3a6", "#accbe6", "#b23000", "#4b008c", "#ced936", 
                            "#307cbf", "#6d00cc", "#33241a", "#7453a6", "#778000", "#ff80d5", 
                            "#8c8a69", "#ff8c40", "#294d00", "#234010", "#8c40ff", "#69238c", 
                            "#7ca3a6", "#4c4359", "#0000bf", "#a6538a", "grey"), 
                            pt.size = 1.0, group.by = "TFP_Annotations", label = FALSE, order = FALSE) + NoLegend() + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"))
ggsave("annotations.png", plot = annotations_plot, width = 18, height = 20, units = "cm")
