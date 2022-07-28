#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.



#loaded the objects FCAloomToSeurat2.rds and FCA_cds_germline_clustered.Robj
#2022_07: The testis FCA SEURAT object can now be created by reading in an uopdated .rds: FCAloomToSeurat2TFP_Annotations.rds

#loaded the following libraries
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(knitr)
library(methods)
library(monocle3)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(Seurat)
library(tidyr)
library(tibble)

#Create Subsets for germline and spermatid clusters
FCA_germline <- subset(FCA_data, idents = c(25, 22, 5, 78, 40, 97, 41, 51, 73, 35, 33, 48, 105, 64, 45, 13, 21, 56, 100, 3, 106, 16, 55, 63, 39, 27, 83, 15, 0, 28, 26, 24, 18, 50, 17, 96, 12, 30, 37, 75, 66, 10, 46))
spermatids <- subset(FCA_data, idents = c(27, 83, 15, 0, 28, 26, 24, 18, 50, 17, 96, 12, 30, 37, 75, 66, 10, 46))

#For ceation of the dot plot in Figure 5G
levels(FCA_germline) <- c(25, 22, 5, 78, 40, 97, 41, 51, 73, 35, 33, 48, 105, 64, 45, 13, 21, 56, 100, 3, 106, 16, 55, 63, 39, 27, 83, 15, 0, 28, 26, 24, 18, 50, 17, 96, 12, 30, 37, 75, 66, 10, 46)
DotPlot(FCA_germline, cols = c("lightgrey","navyblue"), features = c("Glut3","Parp16","CG6701","Pp2C1")) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_size(range = c(0.1, 10))

#For ceation of the expression over pseudotime plots in Figure 5H
spermatid_genes_over_pt <- c("Pp2C1","CG6701","Parp16","Glut3")
spermatid_cds <- FCA_cds_subset[rowData(FCA_cds_subset)$gene_short_name %in% spermatid_genes_over_pt]
plot_genes_in_pseudotime(spermatid_cds, panel_order = c("Pp2C1","CG6701","Parp16","Glut3"), cell_size = 0.005)
#For adding the fold changes
#determine overall max and repeat for each gene
test_cds <- FCA_cds_subset[rowData(FCA_cds_subset)$gene_short_name %in% "Pp2C1"]
test_tid <- plot_genes_in_pseudotime(test_cds, cell_size = 0.005)
max(test_tid$data$expectation)
#make new cds with cutoff around 25 PT
cds_no_spermatids <- choose_cells(FCA_cds_subset)
#repeat to determine the <25 PT max expression (again repeat with each gene)
test_pre_cds <- cds_no_spermatids[rowData(cds_no_spermatids)$gene_short_name %in% "Glut3"]
test_pre_tid <- plot_genes_in_pseudotime(test_pre_cds, cell_size = 0.005)
max(test_pre_tid$data$expectation)
#Fold change was calculated by dividing the overall max by the <25PT max

#Code used for the supplemental figure associated with figure 5
#For creation of the UMAP plots
FeaturePlot(FCA_germline, features = c("soti"), cols = c("lightgrey","navyblue"))
FeaturePlot(FCA_germline, features = c("wa-cup"), cols = c("lightgrey","navyblue"))
FeaturePlot(FCA_germline, features = c("loqs"), cols = c("lightgrey","navyblue"))
FeaturePlot(FCA_germline, features = c("f-cup"), cols = c("lightgrey","navyblue"))

#For the scatter plot
#To generage a table with Pearson CC values for how coorelation each of the 162 spermatid transcribed genes is with each other
PearsonCC_spermatids <- list()
scatterplots_spermatids <- list()
Gene_1 <- c("ACXA","Adck5","Adgf-B","bbc","boly","c-cup","cdc14","CG10512","CG10749","CG10845","CG10969","CG11226","CG11286","CG11298","CG1137","CG11629","CG11779","CG12035","CG12126","CG12209","CG12448","CG12783","CG12983","CG13127","CG13330","CG13337","CG13494","CG13526","CG13538","CG13989","CG14011","CG14106","CG14294","CG14505","CG14546","CG14605","CG14689","CG14739","CG14864","CG15025","CG15059","CG15425","CG15548","CG15638","CG1571","CG15767","CG16825","CG16888","CG16984","CG17098","CG17564","CG17666","CG18266","CG18368","CG1999","CG30099","CG30110","CG30178","CG30284","CG30325","CG30334","CG30392","CG30398","CG30416","CG30484","CG30485","CG31128","CG31204","CG31275","CG31404","CG31730","CG31735","CG31784","CG31797","CG31815","CG31817","CG31897","CG32086","CG32106","CG32117","CG32161","CG32228","CG32396","CG32681","CG32988","CG33017","CG33061","CG33125","CG3408","CG34167","CG34216","CG34432","CG3748","CG4073","CG42266","CG42570","CG42758","CG43339","CG43449","CG43691","CG43796","CG44433","CG44774","CG46385","CG4714","CG5478","CG5561","CG5906","CG6059","CG6136","CG6333","CG6652","CG6701","CG8043","CG8136","CG8654","CG8851","d-cup","Efhc1.2","f-cup","glob3","Glut3","h-cup","hale","hpRNA:CR18854","hpRNA:CR32207","Imp","lncRNA:CR31781","lncRNA:CR43939","lncRNA:CR44225","loqs","m-cup","MCU","mex1","Mst36Fb","Mtl","orb","p-cup","Parp16","PCB","Pglym87","Pif2","Pka-C2","Pp2C1","PpD6","prage","r-cup","Rab3-GAP","schuy","SdhAL","soti","ste24c","sunz","t-cup","TER94","Tob","TTLL1B","Vmat","w-cup","wa-cup","whip","yps")
Gene_2 <- c("ACXA","Adck5","Adgf-B","bbc","boly","c-cup","cdc14","CG10512","CG10749","CG10845","CG10969","CG11226","CG11286","CG11298","CG1137","CG11629","CG11779","CG12035","CG12126","CG12209","CG12448","CG12783","CG12983","CG13127","CG13330","CG13337","CG13494","CG13526","CG13538","CG13989","CG14011","CG14106","CG14294","CG14505","CG14546","CG14605","CG14689","CG14739","CG14864","CG15025","CG15059","CG15425","CG15548","CG15638","CG1571","CG15767","CG16825","CG16888","CG16984","CG17098","CG17564","CG17666","CG18266","CG18368","CG1999","CG30099","CG30110","CG30178","CG30284","CG30325","CG30334","CG30392","CG30398","CG30416","CG30484","CG30485","CG31128","CG31204","CG31275","CG31404","CG31730","CG31735","CG31784","CG31797","CG31815","CG31817","CG31897","CG32086","CG32106","CG32117","CG32161","CG32228","CG32396","CG32681","CG32988","CG33017","CG33061","CG33125","CG3408","CG34167","CG34216","CG34432","CG3748","CG4073","CG42266","CG42570","CG42758","CG43339","CG43449","CG43691","CG43796","CG44433","CG44774","CG46385","CG4714","CG5478","CG5561","CG5906","CG6059","CG6136","CG6333","CG6652","CG6701","CG8043","CG8136","CG8654","CG8851","d-cup","Efhc1.2","f-cup","glob3","Glut3","h-cup","hale","hpRNA:CR18854","hpRNA:CR32207","Imp","lncRNA:CR31781","lncRNA:CR43939","lncRNA:CR44225","loqs","m-cup","MCU","mex1","Mst36Fb","Mtl","orb","p-cup","Parp16","PCB","Pglym87","Pif2","Pka-C2","Pp2C1","PpD6","prage","r-cup","Rab3-GAP","schuy","SdhAL","soti","ste24c","sunz","t-cup","TER94","Tob","TTLL1B","Vmat","w-cup","wa-cup","whip","yps")
gene_pairs <- crossing(Gene_1, Gene_2)
x <- (1:26244)
for (i in x) {
  feature1_gene <- gene_pairs[[i, 1]]
  feature2_gene <- gene_pairs[[i, 2]]
  scatterplots_spermatids[[i]]<-FeatureScatter(spermatids, feature1 = feature1_gene, feature2 = feature2_gene)
  PearsonCC_spermatids[[i]] <- scatterplots_spermatids[[i]]$labels$title
}
PearsonCC_spermatids <- as.data.frame(PearsonCC_spermatids)
PearsonCC_spermatids <- t(PearsonCC_spermatids)
PearsonCC_df_spermatids <- as.data.frame(matrix(PearsonCC_spermatids[,1], ncol=162))
x <- gene_pairs[1:162, 2]
x <- x$Gene_2
rownames(PearsonCC_df_spermatids) <- x
colnames(PearsonCC_df_spermatids) <- x
PearsonCC_df_spermatids[PearsonCC_df_spermatids == 1] <- NA
gene_order <- c("soti","Mst36Fb","CG12209","CG46385","CG34432","CG4714","CG42266","CG43691","CG33017","h-cup","CG30178","CG31815","hale","CG13127","CG32228","CG43339","PpD6","CG15767","CG30325","CG15425","CG30485","CG13526","CG34167","Rab3-GAP","CG12983","CG31784","CG14864","CG30392","CG30398","boly","CG10845","CG1137","CG11779","CG32086","CG17564","CG42570","CG12783","CG13538","CG12035","CG30484","CG31897","CG30334","Pif2","CG31404","CG8851","f-cup","Glut3","CG13330","CG14294","CG18368","CG31735","CG6652","Tob","whip","CG30110","wa-cup","CG31797","CG6333","CG43796","orb","loqs","sunz","c-cup","CG31128","CG4073","CG42758","CG44433","Pp2C1","CG10969","CG14546","CG31817","CG32106","CG12126","Imp","CG14106","d-cup","CG15059","Parp16","CG43449","CG30099","r-cup","lncRNA:CR43939","p-cup","CG15548","CG6701","CG10512","Mtl","Adgf-B","lncRNA:CR44225","CG30284","SdhAL","CG32988","TER94","CG3748","CG14505","CG5478","ste24c","CG8654","glob3","CG32117","prage","CG5906","CG8136","CG13989","CG31204","CG8043","CG32161","CG11298","CG31730","ACXA","w-cup","CG16825","mex1","CG33061","CG6136","hpRNA:CR18854","Vmat","CG44774","CG32681","hpRNA:CR32207","Adck5","CG14689","CG14011","CG15025","CG6059","Efhc1.2","lncRNA:CR31781","CG13494","CG5561","CG16984","CG34216","CG17666","MCU","CG10749","TTLL1B","CG15638","CG16888","cdc14","Pglym87","CG30416","schuy","CG33125","Pka-C2","CG11286","CG14739","t-cup","CG32396","CG11629","m-cup","CG18266","CG13337","CG11226","CG3408","CG12448","CG1571","yps","bbc","PCB","CG31275","CG17098","CG14605","CG1999")
PearsonCC_df_spermatids <- PearsonCC_df_spermatids[, gene_order]
PearsonCC_df_spermatids <- PearsonCC_df_spermatids[gene_order ,]
spermatid_heatmap <- pheatmap(PearsonCC_df_spermatids, cluster_rows = FALSE, cluster_cols = FALSE)

PearsonCC_df_spermatids<-as.data.frame(PearsonCC_df_spermatids)
#Exported as .csv. For each gene, the Pearson CC values were averaged and used in the scatter plot

#To calculate the average expression of each of the 162 spermatid transcribed genes in spermatid clusters
SpermAvg<-AverageExpression(spermatids, features = c("ACXA","Adck5","Adgf-B","bbc","boly","c-cup","cdc14","CG10512","CG10749","CG10845","CG10969","CG11226","CG11286","CG11298","CG1137","CG11629","CG11779","CG12035","CG12126","CG12209","CG12448","CG12783","CG12983","CG13127","CG13330","CG13337","CG13494","CG13526","CG13538","CG13989","CG14011","CG14106","CG14294","CG14505","CG14546","CG14605","CG14689","CG14739","CG14864","CG15025","CG15059","CG15425","CG15548","CG15638","CG1571","CG15767","CG16825","CG16888","CG16984","CG17098","CG17564","CG17666","CG18266","CG18368","CG1999","CG30099","CG30110","CG30178","CG30284","CG30325","CG30334","CG30392","CG30398","CG30416","CG30484","CG30485","CG31128","CG31204","CG31275","CG31404","CG31730","CG31735","CG31784","CG31797","CG31815","CG31817","CG31897","CG32086","CG32106","CG32117","CG32161","CG32228","CG32396","CG32681","CG32988","CG33017","CG33061","CG33125","CG3408","CG34167","CG34216","CG34432","CG3748","CG4073","CG42266","CG42570","CG42758","CG43339","CG43449","CG43691","CG43796","CG44433","CG44774","CG46385","CG4714","CG5478","CG5561","CG5906","CG6059","CG6136","CG6333","CG6652","CG6701","CG8043","CG8136","CG8654","CG8851","d-cup","Efhc1.2","f-cup","glob3","Glut3","h-cup","hale","hpRNA:CR18854","hpRNA:CR32207","Imp","lncRNA:CR31781","lncRNA:CR43939","lncRNA:CR44225","loqs","m-cup","MCU","mex1","Mst36Fb","Mtl","orb","p-cup","Parp16","PCB","Pglym87","Pif2","Pka-C2","Pp2C1","PpD6","prage","r-cup","Rab3-GAP","schuy","SdhAL","soti","ste24c","sunz","t-cup","TER94","Tob","TTLL1B","Vmat","w-cup","wa-cup","whip","yps"))
SpermAvg<-as.data.frame(SpermAvg)
#Exported as .csv and values were averaged for use in the scatter plot
