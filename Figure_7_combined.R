#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.

#  2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.  

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.


# ********* Figure_7_AtoE_code.R **************

#load packages
library(Seurat)
library(tidyverse)
library(SCopeLoomR)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)

memory.limit(size=56000) ##increases the memory limit on a pc
FCA_bioHub_ScrapeD <- readRDS(file = "C:/Users/jasgr/Downloads/FCAloomToSeurat2.rds") #Input Seurat Object
#2022_07: This oibject can now be created by reading in an updated testis FCA SEURAT object: FCAloomToSeurat2TFP_Annotations.rds


#subset out hub ####
Hub <- subset(FCA_bioHub_ScrapeD, idents = c("90"))

#examine genes of interest ###
Hub_genes2 <- c("upd1", "hh", "CadN", "Fas3", "tup") #features for the featureplot

#Making the plots one at a time

##upd1 Figure 7 A
upd1_plot <- FeaturePlot(Hub, features = c("upd1"), pt.size = 1.5, cols = c("lightgrey", "navyblue"), combine = FALSE)
upd1Patch_Plot <- patchwork::wrap_plots(upd1_plot) * ylim(c(9, 10.2)) * xlim(c(6.25, 7.25))

upd1_notitle <- upd1Patch_Plot + labs(title = NULL, x = NULL, y = NULL)
White_upd1 <- upd1_notitle + theme(panel.background = element_rect(fill = "white"),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   panel.border = element_rect(color = "red", size = 2),
                                   axis.line = element_blank())
White_upd1
upd1_Full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("upd1"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
upd1_Whole_Plot <-  upd1_Full_Plot + theme(plot.title = element_text(size = 20))

upd1_TogetherPlot <- upd1_Whole_Plot +
  inset_element(p = White_upd1,
                left = 0.02,
                bottom = 0.7,
                right = 0.3,
                top = 1)
upd1_TogetherPlot
ggsave("upd1_w_inset.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

##hh Figure 7 B
hh_plot <- FeaturePlot(Hub, features = c("hh"), pt.size = 1.5, cols = c("lightgrey", "navyblue"), combine = FALSE)
hhPatch_Plot <- patchwork::wrap_plots(hh_plot) * ylim(c(9, 10.2)) * xlim(c(6.25, 7.25))

hh_notitle <- hhPatch_Plot + labs(title = NULL, x = NULL, y = NULL)
White_hh <- hh_notitle + theme(panel.background = element_rect(fill = "white"),
                               axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               panel.border = element_rect(color = "red", size = 2),
                               axis.line = element_blank())
White_hh
hh_Full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("hh"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
hh_Whole_Plot <-  hh_Full_Plot + theme(plot.title = element_text(size = 20))

hh_TogetherPlot <- hh_Whole_Plot +
  inset_element(p = White_hh,
                left = 0.02,
                bottom = 0.7,
                right = 0.3,
                top = 1)
hh_TogetherPlot
ggsave("hh_w_inset.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

##CadN Figure 7 C
CadN_plot <- FeaturePlot(Hub, features = c("CadN"), pt.size = 1.5, cols = c("lightgrey", "navyblue"), combine = FALSE)
CadNPatch_Plot <- patchwork::wrap_plots(CadN_plot) * ylim(c(9, 10.2)) * xlim(c(6.25, 7.25))

CadN_notitle <- CadNPatch_Plot + labs(title = NULL, x = NULL, y = NULL)
WhiteCadN <- CadN_notitle + theme(panel.background = element_rect(fill = "white"),
                                  panel.border = element_rect(colour = "red", size=2),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.line = element_blank())

WhiteCadN
CadN_full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("CadN"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
CadN_Whole_Plot <- CadN_full_Plot + theme(plot.title = element_text(size = 20))
CadN_Whole_Plot

CadN_TogetherPlot <- CadN_Whole_Plot +
  inset_element(p = WhiteCadN,
                left = 0.02,
                bottom = 0.7,
                right = 0.3,
                top = 1)
CadN_TogetherPlot
ggsave("CadN_w_inset.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

##Fas3 Figure 7 D
Fas3_plot <- FeaturePlot(Hub, features = c("Fas3"), pt.size = 1.5, cols = c("lightgrey", "navyblue"), combine = FALSE)
Fas3Patch_Plot <- patchwork::wrap_plots(Fas3_plot) * ylim(c(9, 10.2)) * xlim(c(6.25, 7.25))

Fas3_notitle <- Fas3Patch_Plot + labs(title = NULL, x = NULL, y = NULL)
WhiteFas3 <- Fas3_notitle + theme(panel.background = element_rect(fill = "white"),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  panel.border = element_rect(color = "red", size = 2),
                                  axis.line = element_blank())
WhiteFas3
Fas3_Full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("Fas3"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
Fas3_Whole_Plot <- Fas3_Full_Plot + theme(plot.title = element_text(size = 20))

Fas3_TogetherPlot <- Fas3_Whole_Plot +
  inset_element(p = WhiteFas3,
                left = 0.02,
                bottom = 0.7,
                right = 0.3,
                top = 1)
Fas3_TogetherPlot
ggsave("Fas3_w_inset.png", plot = last_plot(), width = 7, height = 7, dpi = 300)

##tup Figure 7 E
tup_plot <- FeaturePlot(Hub, features = c("tup"), pt.size = 1.5, cols = c("lightgrey", "navyblue"), combine = FALSE)
tupPatch_Plot <- patchwork::wrap_plots(tup_plot) * ylim(c(9, 10.2)) * xlim(c(6.25, 7.25))

tup_notitle <- tupPatch_Plot + labs(title = NULL, x = NULL, y = NULL)
White_tup <- tup_notitle + theme(panel.background = element_rect(fill = "white"),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 panel.border = element_rect(color = "red", size = 2),
                                 axis.line = element_blank())
White_tup
tup_Full_Plot <- FeaturePlot(FCA_bioHub_ScrapeD, features = c("tup"), pt.size = 1.0, cols = c("lightgrey", "navyblue"))
tup_Whole_Plot <- tup_Full_Plot + theme(plot.title = element_text(size = 20))

tup_TogetherPlot <- tup_Whole_Plot +
  inset_element(p = White_tup,
                left = 0.02,
                bottom = 0.7,
                right = 0.3,
                top = 1)
tup_TogetherPlot
ggsave("tup_w_inset.png", plot = last_plot(), width = 7, height = 7, dpi = 300)


# ******** Figure_7_DE_code.R***********

library(Seurat)
library(SeuratData)
library(BiocManager)
library(limma)

FCA_bioHub_ScrapeD <- readRDS(file = "C:/Users/jasgr/Downloads/FCAloomToSeurat2.rds")
#2022_07: This oibject can now be created by reading in an updated testis FCA SEURAT object: FCAloomToSeurat2TFP_Annotations.rds

levels(FCA_bioHub_ScrapeD)

#Hub markers
Hub <- FindMarkers(FCA_bioHub_ScrapeD, ident.1 = "90", ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1.0)
head(Hub)
write.table(Hub, file = "Hub_markers.csv", sep=',')

#CySCs markers
CySCs_markers <- FindMarkers(FCA_bioHub_ScrapeD, ident.1 = "62", ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1.0)
head(CySCs_markers)
write.table(CySCs_markers, file = "CySCs_markers.csv", sep=',')

#Gonia markers
Gonia_markers <- FindMarkers(FCA_bioHub_ScrapeD, ident.1 = "25", ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1.0)
head(Gonia_markers)
write.table(Gonia_markers, file = "Gonia_markers.csv", sep=',')

#Terminal Epithelium Markers
TE_markers <- FindMarkers(FCA_bioHub_ScrapeD, ident.1 = c("42", "20", "89"), ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1.0)
head(TE_markers)
write.table(TE_markers, file = "TE_markers.csv", sep=',')

#Hemocytes markers
Hemocytes_markers <- FindMarkers(FCA_bioHub_ScrapeD, ident.1 = c("44", "52", "109"), ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1.0)
head(Hemocytes_markers)
write.table(Hemocytes_markers, file = "Hemocytes_markers.csv", sep=',')


# ********** Figure_7_VennDiagram_code.R ***********

library("gdata")

#Followed this protocol: https://www.r-bloggers.com/2016/06/working-with-venn-diagrams/
#import list of up-regulated hub genes and up-regulated CySC genes and look at the head and tail
head(HubandCySCs_deup1)
tail(HubandCySCs_deup1)

#To convert the data frame to separate gene lists with the empty strings removed use lapply() with the function(x) x[x != ""]
HubCyDE <- lapply(as.list(HubandCySCs_deup1), function(x) x[x != ""])
lapply(HubCyDE, tail)

#Make Venn Diagram
require("VennDiagram")
venn.plot.hubCy <- venn.diagram(HubCyDE, NULL, cex = 2, main="Hub vs CySCs")

#Need to use grid.newpage() before every Venn Diagram or they will just stack on top of each other
grid.newpage()
grid.draw(venn.plot.hubCy)


#import list of up-regulated hub genes and up-regulated Gonia genes
head(HubandGonia_deup1)
tail(HubandGonia_deup1)

#To convert the data frame to separate gene lists with the empty strings removed use lapply() with the function(x) x[x != ""]
HubSGDE <- lapply(as.list(HubandGonia_deup1), function(x) x[x != ""])
lapply(HubSGDE, tail)

#Make Venn Diagram
require("VennDiagram")
venn.plot.hubSG <- venn.diagram(HubSGDE, NULL, cex = 2, main="Hub vs Gonia")

#Need to use grid.newpage() before every Venn Diagram or they will just stack on top of each other
grid.newpage()
grid.draw(venn.plot.hubSG)


#import list of up-regulated hub genes and up-regulated Hemocyte genes
head(HubandHemocytes_deup1)
tail(HubandHemocytes_deup1)

#To convert the data frame to separate gene lists with the empty strings removed use lapply() with the function(x) x[x != ""]
HubH_DE <- lapply(as.list(HubandHemocytes_deup1), function(x) x[x != ""])
lapply(HubH_DE, tail)

#Make Venn Diagram
require("VennDiagram")
venn.plot.hubH <- venn.diagram(HubH_DE, NULL, cex = 2, main="Hub vs Hemocytes")

#Need to use grid.newpage() before every Venn Diagram or they will just stack on top of each other
grid.newpage()
grid.draw(venn.plot.hubH)


#import list of up-regulated hub genes and up-regulated Terminal Epithelia genes
head(HubandTE_deup1)
tail(HubandTE_deup1)

#To convert the data frame to separate gene lists with the empty strings removed use lapply() with the function(x) x[x != ""]
HubTE_DE <- lapply(as.list(HubandTE_deup1), function(x) x[x != ""])
lapply(HubTE_DE, tail)

#Make Venn Diagram
require("VennDiagram")
venn.plot.hubTE <- venn.diagram(HubTE_DE, NULL, cex = 2, main="Hub vs TE")

#Need to use grid.newpage() before every Venn Diagram or they will just stack on top of each other
grid.newpage()
grid.draw(venn.plot.hubTE)

##The following is code to get the list of genes in each segment of the Venn Diagram

require("gplots")
#Hub vs CySCs
HubCy_together <- venn(HubCyDE, show.plot=FALSE)
# You can inspect the contents of this object with the str() function
str(HubCy_together)

intersCy <- attr(HubCy_together,"intersections")
lapply(intersCy, head) 

write.table(intersCy$CySCs, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/CySCs_Venn_Files/CySCsvsHub_CySCsonly.csv", sep=',')
write.table(intersCy$Hub, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/CySCs_Venn_Files/CySCsvsHub_Hubonly.csv", sep=',')
write.table(intersCy$`Hub:CySCs`, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/CySCs_Venn_Files/CySCsvsHub_sharedgenes.csv", sep=',')

#Hub vs Gonia
HubSG_together <- venn(HubSGDE, show.plot=FALSE)
# You can inspect the contents of this object with the str() function
str(HubSG_together)
intersSG <- attr(HubSG_together,"intersections")
lapply(intersSG, head) 

write.table(intersSG$Gonia, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/Gonia_Venn_Files/GoniavsHub_Goniaonly.csv", sep=',')
write.table(intersSG$Hub, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/Gonia_Venn_Files/GoniavsHub_Hubonly.csv", sep=',')
write.table(intersSG$`Hub:Gonia`, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/Gonia_Venn_Files/GoniavsHub_sharedgenes.csv", sep=',')

#Hub vs Hemocytes
HubHemo_together <- venn(HubH_DE, show.plot=FALSE)
# You can inspect the contents of this object with the str() function
str(HubHemo_together)
intersHemo <- attr(HubHemo_together,"intersections")
lapply(intersHemo, head) 

write.table(intersHemo$Hemocytes, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/Hemocytes_Venn_Files/HemovsHub_Hemoonly.csv", sep=',')
write.table(intersHemo$Hub, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/Hemocytes_Venn_Files/HemovsHub_Hubonly.csv", sep=',')
write.table(intersHemo$`Hub:Hemocytes`, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/Hemocytes_Venn_Files/HemovsHub_sharedgenes.csv", sep=',')


# Hub vs TE
HubTE_together <- venn(HubTE_DE, show.plot=FALSE)
# You can inspect the contents of this object with the str() function
str(HubTE_together)

intersTEHub <- attr(HubTE_together,"intersections")
lapply(intersTEHub, head) 
write.table(intersTEHub$TE, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/TE_Venn_Files/TEvsHub_TEonlygenes.csv", sep=',')
write.table(intersTEHub$`Hub:TE`, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/TE_Venn_Files/TEvsHub_sharedgenes.csv", sep=',')
write.table(intersTEHub$Hub, file = "N:/Jasmine Grey/Fly Cell Atlas/Final_R_TFP_files/TFP_R_files_final/TE_Venn_Files/TEvsHub_Hubonlygenes.csv", sep=',')
