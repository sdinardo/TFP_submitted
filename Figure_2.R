#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.

#  2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.  

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.



##load packages ####
library(Seurat)
library(tidyverse)
library(SCopeLoomR)
library(Matrix)
library(ggplot2)
library(dplyr)

setwd("~/testisseurat")

memory.limit(size=56000) ##increases the memory limit on a pc
FCA_bioHub_ScrapeB <- readRDS(file = "/Users/mara/Documents/TestisSeurat/FCA_bioHub_ScrapeB.rds")

#subset out gonia ####
gonia3 <- subset(FCA_bioHub_ScrapeB, idents = c("25", "22", "5", "78", "40", "97", "41", "51", "73", "35", "33", "48", "105", "21", "56", "100", "3", "106", "16", "64", "45", "13", "55", "63", "39", "27", "83", "15", "0", "28", "26", "24", "18", "50", "17", "96", "12", "30", "37", "75", "66", "10", "46"))
#examine genes of interest ####
gl_progression_features <- c( "stg","zpg", "bam", "aub", "kmg",  "Rbp4", "aly", "can", "sa", "CycB", "fzo", "Dic61B", "p-cup") #features for the dotplot

gl_progression_features_4 <- c("aub", "kmg", "Rbp4", "zpg", "CycB", "fzo", "p-cup", "aly") #features for the featureplot


# save each individual feature plot for each gene in gl_progression_features_4
for (gene in gl_progression_features_4) {
  FeaturePlot(gonia3, features = c(gene), pt.size = 0.5, cols = c("lightgrey", "navyblue"))
  ggsave(filename = paste(gene, ".jpg"), plot=last_plot(), height=4, width=4.5, dpi=600)
  print(gene)
}


# Visualize co-expression of two features simultaneously
FeaturePlot(gonia3, features = c("aub", "fzo"), blend = TRUE, cols = c("blue", "red"), pt.size = 1)
ggsave("aub fzo.jpg",  plot=last_plot(), height=3.2, width=12, dpi=600)


# To set the order of clusters from gonia --> SpCyte:
levels(gonia3) <- c(25, 22, 5, 78, 40, 97, 41, 51, 73, 35, 33, 48, 105, 64, 45, 13, 21, 56, 100, 3, 106, 16, 55, 63, 39, 27, 83, 15, 0, 28, 26, 24, 18, 50, 17, 96, 12, 30, 37, 75, 66, 10, 46)

# Return the DotPlot ####
DotPlot(gonia3, features = rev(gl_progression_features), cols = c("lightgrey", "navyblue")) + RotatedAxis()

# Annotate this plot in ggplpot
glProgressionDotPlot <- DotPlot(gonia3, features = rev(gl_progression_features)) + RotatedAxis()

glProgressionDotPlot + labs(title = "Expression by cluster as germline progresses") +
  xlab("Germline Markers") +
  ylab("") + 
  theme(plot.title=element_text(size=15)) + coord_flip() +  ##coord_flip() swaps the axes and rev() reverses the order of a vector
  scale_color_gradient(low = "lightgrey", high = "navyblue") +
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13))

spermatogonia<-c("25","22")
spermatocyte<-c("5","40","41","51","35","33","48","21","3","13","16","39","45","55","56","64","73","78","97","100","105","106")
spermatid<-c("30","0","15","18","24","26","28","12","17","27","50","63","83","96","10","37","46","66","73")

#save dotplot
ggsave("gonia dotplot with bam.jpg", plot=last_plot(), height = 4, width = 14, dpi=600)

#make UMAP with cluster labels
DimPlot(gonia3, label = TRUE, label.size = 4) + NoLegend()

ggsave("germline umap.jpg", plot = last_plot(), scale = 1, dpi = 600)
