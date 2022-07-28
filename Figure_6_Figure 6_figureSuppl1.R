#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.

#  2022_07_28: The current R data file, FCAloomToSeurat2TFP_Annotations.rds, can be read in to produce the Seurat Object used in Raz et al.  

# Users may need to edit code to their purpose, e.g., changing filePaths and fileNames as necessary.



#The following codes were used to make figure 6 (Progression of differentiation in the somatic cyst cell lineage) for the Testis follow-up paper to the Fly Cell Atlas paper
#Use testis FCA SEURAT object readRDS 'FCAloomToSeurat2.rds'
  # 1.This object has all the info from the loom
  # 2.The counts have been scaled and normalized
  # 3.The definitive hub cluster nuclei have been correctly allocated.
#2022_07: This oibject can now be created by reading in an uopdated testis FCA SEURAT object: FCAloomToSeurat2TFP_Annotations.rds

# load packages   ####
library(Seurat)
library(tidyverse)
library(SCopeLoomR)
library(Matrix)
library(plotly)

FCA_loomToSeurat <- readRDS(file = "~/Desktop/FCA R/FCAloomToSeurat2.rds")

# Check UMAP of entire testis
DimPlot(FCA_loomToSeurat, label = TRUE) + NoLegend()  # Clusters Leiden 6.0 with wings separated from hub cluster

#Subset out cyst lineage
cyst <- subset(FCA_loomToSeurat, idents = c("62", "36", "58", "77", "65", "47","95","88","57","80","84","72","74","98","69","67","60","94","104","81","79","68"))

## Create UMAP of cyst lineage for figure 6A##
my_cols_cyst <- c("62"="#B185FE","36"="#E57F36", "58"="#ACA100", "77"="#62B302","65"="#00BD5B","47"="#8EAA00","88"="#0CC0A7","57"="#DA8E00","80"="#0BBEC5","84"="#7E95FF","72"="#00A5FF","74"="#06BADF","98"="#FF6996","69"="#00C285","67"="#F8766D","60"="#00B2F3","94"="#EF67EC","104"="#D972FF","81"="#FD61D3","79"="#CB9501","68"="#FF64B6","95"="#E0EF27")
DimPlot(cyst, label = TRUE, label.size = 4, pt.size = 0.2, cols=my_cols_cyst) + NoLegend() 

## Visualize co-expression of two features, tj, eya, and so, simultaneously for figure 6B. tj&eya and eya&so maps merged in imageJ by making stack --> Z project --> min intensity##
## Separated FeaturePlots included in supplemental material (Figure 6 - supplement 1 C, D) with corresponding heatmaps##
FeaturePlot(cyst, features = c("tj", "eya"), blend = TRUE, cols = c("blue", "orange"), pt.size = 1, combine = FALSE)
tjeya <- FeaturePlot(cyst, features = c("tj", "eya"), blend = TRUE, cols = c("blue", "orange"), pt.size = 1, combine = FALSE)
tjeya[[3]] + FontSize(x.title = 12, y.title = 12) + NoLegend() + labs(title = "")
tjeya[[4]] + FontSize(x.title = 20, y.title = 20) + labs(title = "")

FeaturePlot(cyst, features = c("eya", "so"), blend = TRUE, cols = c("orange", "red"), pt.size = 1, combine = FALSE)
eyaso <- FeaturePlot(cyst, features = c("eya", "so"), blend = TRUE, cols = c("orange", "red"), pt.size = 1, combine = FALSE)
eyaso[[3]] + FontSize(x.title = 12, y.title = 12) + NoLegend() + labs(title = "")
eyaso[[4]] + FontSize(x.title = 20, y.title = 20) + labs(title = "")

## Visualize expression of individual features for figures 6F, J, M, P, R, T##
cl_progression_features <- c("piwi", "Amph", "geko", "eya", "shg", "Nep4", "Akr1B", "CG3902")
for (gene in cl_progression_features) {
  FeaturePlot(cyst, features = c(gene), pt.size = 0.5, cols = c("lightgrey","navyblue")) 
  ggsave(filename = paste(gene, ".jpg"), plot=last_plot(), height=4, width=4.5, dpi=600)
  print(gene)
}

# To set the order of clusters from early --> late Cyst lineage
levels(cyst) <- c("62", "36", "58", "81", "77", "65", "98", "47", "104", "88", "94", "95", "74", "67", "60", "69", "79", "68", "57", "80", "84", "72")

## Analyze the Average UMI per cluster during cyst lineage progression for figure 6E##
genes_UMI_byCluster_cl <- cyst@meta.data %>%
  group_by(leiden_res_6.0, nCount_RNA, nFeature_RNA) # Extract 3 columns meta.data

cl_testPlot <- genes_UMI_byCluster_cl %>%
  group_by(leiden_res_6.0) %>%
  summarise(aveUMI = mean(nCount_RNA), aveGenes = mean(nFeature_RNA)) # Calculate means

# Make the cluster column a character vector
cl_testPlot$leiden_res_6.0 <- as.character(cl_testPlot$leiden_res_6.0)
# Then turn it back into a factor with the levels' in the correct order
cl_testPlot$leiden_res_6.0 <- factor(cl_testPlot$leiden_res_6.0, levels=c("62", "36", "58", "81", "77", "65", "98", "47", "104", "88", "94", "95", "74", "67", "60", "69", "79", "68", "57", "80", "84", "72"))

# Create plot with labels
ggplot(cl_testPlot, aes(leiden_res_6.0,aveUMI)) +
  geom_point(aes(size=aveGenes)) +
  labs(title = "Average UMI per cluster during cyst lineage progression") +
  xlab("CySC to late cyst lineage") +
  ylab("Avearage UMI") + 
  theme(plot.title = element_text(size = 15)) + 
  FontSize(x.title = 14, y.title = 15) 

## Create a DotPlot for the cyst lineage with key features for figure 6V##
cl2_progression_features <- c("tj", "CG3902", "piwi", "Amph", "eya", "Akr1B", "geko", "shg", "so", "Nep4")

DotPlot(cyst, features = rev(cl2_progression_features), cols = c("lightgrey", "navyblue")) + RotatedAxis() 

# Annotate this plot in ggplpot
clProgressionDotPlot <- DotPlot(cyst, features = rev(cl2_progression_features), cols = c("lightgrey", "navyblue")) + RotatedAxis()

clProgressionDotPlot + labs(title = "Expression by cluster as cyst lineage progresses") +
xlab("Cyst Lineage Markers") + ylab("CySC to late cyst lineage") + coord_flip() + theme(text=element_text(size = 9)) + theme(plot.title = element_text(size = 16)) + scale_color_gradient(low = "lightgrey", high = "navyblue") + scale_size(range = c(1, 7)) + FontSize(x.title = 15, y.title = 15)
