#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.

#  2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.  

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.


#load packages
library(Seurat)
library(tidyverse)
library(SCopeLoomR)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)

memory.limit(size=56000) ##increases the memory limit on a pc
FCA_bioHub_ScrapeD <- readRDS(file = "C:/Users/jasgr/Downloads/FCAloomToSeurat2.rds")



#Figure 7 - figure supplement 1A
Original_Clust90 <- DimPlot(FCA_bioHub_ScrapeD, label=F, cells.highlight=WhichCells(FCA_bioHub_ScrapeD, idents = c("90", "111")), pt.size = 1.0, cols.highlight = c("navyblue")) + NoLegend()
Original_Clust90
ggsave("Org_clust90.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

#Figure 7 - figure supplement 1B
rbp4_full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("Rbp4"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
rbp4_Whole_Plot <- rbp4_full_Plot + theme(plot.title = element_text(size = 20))
rbp4_Whole_Plot
ggsave("Rbp4_plot.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

#Figure 7 - figure supplement 1C
zpg_full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("zpg"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
zpg_Whole_Plot <- zpg_full_Plot + theme(plot.title = element_text(size = 20))
zpg_Whole_Plot
ggsave("zpg_plot.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

#Figure 7 - figure supplement 1D
p53_full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("p53"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
p53_Whole_Plot <- p53_full_Plot + theme(plot.title = element_text(size = 20))
p53_Whole_Plot
ggsave("p53_plot.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

#Figure 7 - figure supplement 1E
vas_full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("vas"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
vas_Whole_Plot <- vas_full_Plot + theme(plot.title = element_text(size = 20))
vas_Whole_Plot
ggsave("vas_plot.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

