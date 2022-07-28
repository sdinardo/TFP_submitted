#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.

#  2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.  

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.


#The following codes were used to make Figure 8 (Epithelial Cells of the Testis Organ) and Figure 8 - figure supplement 1 for the Testis follow-up paper to the Fly Cell Atlas paper


# load packages   ####
library(Seurat)
library(tidyverse)
library(SCopeLoomR)
library(Matrix)

FCA_loomToSeurat <- readRDS(file = "~/Desktop/FCA R/FCAloomToSeurat2.rds")

# Check DimPlot
DimPlot(FCA_loomToSeurat, label = TRUE) + NoLegend()  # Clusters Leiden 6.0 with wings separated from hub cluster

#Subset out basal epithelia
basalepithelia <- subset(FCA_loomToSeurat, idents = c("42", "20", "11", "76", "14", "38", "23", "93", "101", "53", "32", "34", "29", "7", "9", "1", "4", "2", "108", "89"))

## Create UMAP of cyst lineage for Figure 8A##
DimPlot(basalepithelia, label = TRUE, label.size = 4, pt.size = 0.2) + NoLegend() 

## Visualize expression of individual features for Figures 8B, D, F##
basalepithelia_progression_features <- c("hh","MtnA","Nep5")
for (gene in basalepithelia_progression_features) {
  FeaturePlot(basalepithelia, features = c(gene), pt.size = 0.5, cols = c("lightgrey", "navyblue"))
  ggsave(filename = paste(gene, ".jpg"), plot=last_plot(), height=4, width=3.5, dpi=600)
  print(gene)
}

## Visualize expression of individual features using the counts slot for Figure 8 - figure supplement 1A,B,C##
for (gene in basalepithelia_progression_features) {
  FeaturePlot(basalepithelia, features = c(gene), slot = "counts", order = TRUE, pt.size = 0.5, cols = c("lightgrey", "navyblue"))
  ggsave(filename = paste(gene, " count2.jpg"), plot=last_plot(), height=4, width=3.5, dpi=600)
  print(gene)
}

## Visualize co-expression of two features, Nep4 and Nep5, simultaneously for Figure 8H##
FeaturePlot(FCA_loomToSeurat, features = c("Nep4", "Nep5"), blend = TRUE, cols = c("blue", "red"), pt.size = 0.5, combine = FALSE)
nep <- FeaturePlot(FCA_loomToSeurat, features = c("Nep4", "Nep5"), blend = TRUE, cols = c("blue", "red"), pt.size = 0.5, combine = FALSE)
nep[[3]] + FontSize(x.title = 12, y.title = 12) + NoLegend() + labs(title = "")
nep[[4]] + FontSize(x.title = 20, y.title = 20) + labs(title = "")

