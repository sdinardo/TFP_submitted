# Load packages and set working directory
library(Seurat)
library(tidyverse)
library(ggplot2)
library(viridis)

setwd("~/Desktop/Lab/FCA/Follow-Up/Code")

# Import Seurat object created by Gabby Vida and Steve DiNardo
testisv1 <- readRDS(file = "FCAloomToSeurat2TFP_Annotations.rds")

#### USING SEURAT OBJECT ####

# FIGURE 3 S1A - Testis
# This code was written by Steve DiNardo and used to produce the
# panels from the other datasets in the remaining three panels

# Average number of UMIs and genes per clusters
# Extract 3 columns meta.data
genes_UMI_byCluster <- testisv1@meta.data %>%
  group_by(Annotation, nCount_RNA, nFeature_RNA)

# Calculate means
testPlot <- genes_UMI_byCluster %>%
  group_by(Annotation) %>%
  summarise(aveUMI = mean(nCount_RNA), aveGenes = mean(nFeature_RNA))

# Make the cluster column a numeric vector
testPlot$Annotation <- as.numeric(as.character(testPlot$Annotation))

# Plot with labels
testis_umi_plot <- ggplot(testPlot, aes(Annotation, aveUMI)) +
  geom_point(aes(size = aveGenes)) +
  labs(title = "Testis: UMI (ave) per cluster", subtitle = "(aveGenes per cluster, see key)") +
  xlab("Cluster #") +
  ylab("UMI per cluster (ave)") +
  ylim(0, 30000) +
  theme(axis.title = element_text(size = 24), plot.title = element_text(size = 30), axis.text = element_text(size = 18)) +
  scale_size(range = c(4, 9))

ggsave("testis_umi_plot.png", plot = testis_umi_plot, width = 40, height = 15, units = "cm")

# FIGURE 3 S1B
# Code was written by Amelie Raz to produce a greyscale FeaturePlot,
# plotting nFeature_RNA using the scRNA-seq dataset.

# FIGURE 3 S1C
# Code was written by Amelie Raz to produce a VlnPlot,
# plotting nCount_RNA by Identity using the scRNA-seq dataset.

# FIGURE 3 S1D
# Note that this code uses an alternate assay of the Seurat object called "log.counts", one containing log-transformed raw counts
# This plot can be created by switching default assay: DefaultAssay(testisv1) <- "log.counts"
# or by using the Key (used below) to override the default "RNA" assay in specifying the feature argument
Syt1_plot <- FeaturePlot(testisv1, features = "logcounts_Syt1", max.cutoff = "q99", pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15)) + scale_color_viridis(option = "plasma")
ggsave("Syt1.png", plot = Syt1_plot, width = 20, height = 20, units = "cm")

# FIGURE 3 S1E
# Note that this code uses an alternate assay of the Seurat object called "log.counts", one containing log-transformed raw counts
# This plot can be created by switching default assay: DefaultAssay(testisv1) <- "log.counts"
# or by using the Key (used below) to override the default "RNA" assay in specifying the feature argument
grh_plot <- FeaturePlot(testisv1, features = "logcounts_grh", pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15)) + scale_color_viridis(option = "plasma")
ggsave("grh.png", plot = grh_plot, width = 20, height = 20, units = "cm")

# FIGURE 3 S1F
# Note that this code uses an alternate assay of the Seurat object called "log.counts", one containing log-transformed raw counts
# This plot can be created by switching default assay: DefaultAssay(testisv1) <- "log.counts"
# or by using the Key (used below) to override the default "RNA" assay in specifying the feature argument
eya_plot <- FeaturePlot(testisv1, features = "logcounts_eya", pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15)) + scale_color_viridis(option = "plasma")
ggsave("eya.png", plot = eya_plot, width = 20, height = 20, units = "cm")

# FIGURE 3 S1G
# Note that this code uses an alternate assay of the Seurat object called "log.counts", one containing log-transformed raw counts
# This plot can be created by switching default assay: DefaultAssay(testisv1) <- "log.counts"
# or by using the Key (used below) to override the default "RNA" assay in specifying the feature argument
upd1_plot <- FeaturePlot(testisv1, features = "logcounts_upd1", pt.size = 1.0, order = TRUE) + theme(text = element_text(size = 30), axis.text = element_text(size = 30), axis.ticks = element_line(size = 1), axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, "cm"), legend.text = element_text(size = 15)) + scale_color_viridis(option = "plasma")
ggsave("upd1.png", plot = upd1_plot, width = 20, height = 20, units = "cm")

# FIGURE 3 S1I
# Code was written by Amelie Raz to produce two DotPlots,
# plotting two genes (per plot) by Identity using the scRNA-seq dataset.
