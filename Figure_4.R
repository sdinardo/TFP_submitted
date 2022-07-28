#  Code to produce, interrogate and/or manipulate the Seurat object used in “Emergent dynamics of adult stem cell lineages from single nucleus and single cell RNA-Seq of Drosophila testes (Raz et al., https://doi.org/10.1101/2022.07.26.501581).  Code provided ‘as is’, under terms of the MIT license.


#load necessary packages

library(ggplot2)
library(monocle3)
library(jsonlite)
library(tidyverse)
library(SCopeLoomR)
library(Matrix)
library(dplyr)
library(Seurat)
library(viridis)


#use the convert function (from SP) to generate cds object from loom file
FCA_cds = convert("/lab/solexa_yamashita/people/Amelie/FCA_data/r_fca_biohub_testis_10x.loom")

#normalize the data by log factor and calculate dimensionality reduction. 50 dimensions determined empirically 
FCA_cds <- preprocess_cds(FCA_cds, num_dim = 50)

#visualize variance
plot_pc_variance_explained(FCA_cds)

#UMAP dimensionality reduction
FCA_cds <- reduce_dimension(FCA_cds)
plot_cells(FCA_cds)

#visualize diagnostic genes and cluster
plot_cells(FCA_cds, genes=c("aub","vas","soti","p-cup"))
FCA_cds = cluster_cells(FCA_cds, resolution=1e-5)
plot_cells(FCA_cds)

#learn trajectory graph
FCA_cds <- learn_graph(FCA_cds)
plot_cells(FCA_cds)

#calculate pseudotime
FCA_cds <- order_cells(FCA_cds)
plot_cells(FCA_cds, color_cells_by = "pseudotime")

#visualize pseudotime parameters in UMAP projection in Seurat (Fig 4B)
pseudotime_look <- FCA_cds@principal_graph_aux$UMAP$pseudotime
FCA_PTanno <- FCA_bioHub_ScrapeB
FCA_PTanno$pseudotime <- pseudotime_look

#subset out germline based on leiden 6 identity
FCA_PTanno_germline <- subset(FCA_PTanno, idents = c("25","22", "5", "78",  "40", "41", "97", "51", "73", "35", "45", "64", "33", "48","13", "56", "105", "106", "21", "55", "16", "39", "100", "3", "63", "30", "66", "27","0","15","28","26","18","17","12","37","75","10","46","50","83","24","96"))
save(FCA_PTanno_germline, file ="/lab/solexa_yamashita/people/Amelie/FCA_data/FCA_PTanno_germline.Robj")

pdf(file = "/lab/solexa_yamashita/people/Amelie/FCA_data/Fig4B.pdf", width = 4.1, height = 4)
FeaturePlot(FCA_PTanno_germline, features = "pseudotime") +   scale_color_viridis(option="plasma")
dev.off()

#get cells with PT values in Seurat object (all germline cells) and apply to monocle objects
FCA_PTann_PT <- as.data.frame(FCA_PTanno_germline$pseudotime)
FCA_PTann_PT <- FCA_PTann_PT %>% arrange(FCA_PTann_PT)
PT_cells <- rownames(FCA_PTann_PT)
PT_cells <- PT_cells[1:21061]
FCA_cds_subset <- FCA_cds[,PT_cells]
save(FCA_cds_subset, file = "/lab/solexa_yamashita/people/Amelie/FCA_data/FCA_cds_subset_filtered.Robj")

#cluster germline object
FCA_cds_subset = cluster_cells(FCA_cds_subset, resolution = 0.0002)
plot_cells(FCA_cds_subset, show_trajectory_graph = FALSE,group_label_size = 7)
colData(FCA_cds_subset)$cluster_id <- as.character(clusters(FCA_cds_subset))

#Apply Monocle clusters to Seurat object
clusterid_look <- FCA_cds_subset@colData$cluster_id
coldata_look <- FCA_cds_subset@colData@rownames
names(clusterid_look) <- coldata_look

A<-names(which(clusterid_look == "8"))
B<-names(which(clusterid_look == "7"))
C<-names(which(clusterid_look == "6"))
D<-names(which(clusterid_look == "10"))
E<-names(which(clusterid_look == "11"))
Fx<-names(which(clusterid_look == "13"))
G<-names(which(clusterid_look == "12"))
H<-names(which(clusterid_look == "4"))
I<-names(which(clusterid_look == "5"))
J<-names(which(clusterid_look == "9"))
K<-names(which(clusterid_look == "2"))
L<-names(which(clusterid_look == "1"))
M<-names(which(clusterid_look == "14"))
N<-names(which(clusterid_look == "3"))

FCA_germline_anno <- FCA_PTanno_germline

FCA_germline_anno <- SetIdent(FCA_germline_anno, N,"N")
FCA_germline_anno <- SetIdent(FCA_germline_anno, M,"M")
FCA_germline_anno <- SetIdent(FCA_germline_anno, L,"L")
FCA_germline_anno <- SetIdent(FCA_germline_anno, K,"K")
FCA_germline_anno <- SetIdent(FCA_germline_anno, J,"J")
FCA_germline_anno <- SetIdent(FCA_germline_anno, I,"I")
FCA_germline_anno <- SetIdent(FCA_germline_anno, H,"H")
FCA_germline_anno <- SetIdent(FCA_germline_anno, G,"G")
FCA_germline_anno <- SetIdent(FCA_germline_anno, Fx,"F")
FCA_germline_anno <- SetIdent(FCA_germline_anno, E,"E")
FCA_germline_anno <- SetIdent(FCA_germline_anno, D,"D")
FCA_germline_anno <- SetIdent(FCA_germline_anno, C,"C")
FCA_germline_anno <- SetIdent(FCA_germline_anno, B,"B")
FCA_germline_anno <- SetIdent(FCA_germline_anno, A,"A")

colors = c('#F8766D', '#DB8E01', '#AFA200','#64B201','#64B201','#02BD5D','#02C1A7','#02C1A7','#02C1A7','#02C1A7','#02BADE','#00A6FF','#EF67EB','#FF62B6','#7F7F7F','#7F7F7F','#7F7F7F','#7F7F7F')
pdf(file = "/lab/solexa_yamashita/people/Amelie/FCA_data/Fig4A.pdf", width = 4.2, height = 4)
DimPlot(FCA_germline_anno, cols = colors)
dev.off()
save(FCA_germline_anno, file = "/lab/solexa_yamashita/people/Amelie/FCA_data/FCA_germline_anno.Robj")


#####
#generate monocle object for the cell data
cds <- load_cellranger_data("/home/araz/R/10x_data")
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
plot_cells(cds)

#generate and display germline pseudotime trajectory
cds <- learn_graph(cds)
plot_cells(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime")
cds_subset <- choose_cells(cds)
plot_cells(cds_subset, color_cells_by = "pseudotime", show_trajectory_graph = FALSE)

#cluster germline
cds_subset = cluster_cells(cds_subset, resolution = 0.003)

#annotate clusters
plot_cells(cds_subset, color_cells_by = "cluster_id", show_trajectory_graph = FALSE)
colData(cds_subset)$cluster_id <- as.character(clusters(cds_subset))
colData(cds_subset)$cluster_id = dplyr::recode(colData(cds_subset)$cluster_id,
                                               "1"="F",
                                               "2"="D",
                                               "3"="F",
                                               "4"="E",
                                               "5"="F",
                                               "6"="B",
                                               "7"="H",
                                               "8"="C",
                                               "9"="H",
                                               "10"="E",
                                               "11"="I",
                                               "12"="G",
                                               "13"="G",
                                               "14"="I",
                                               "15"="I",
                                               "16"="K",
                                               "17"="G",
                                               "18"="H",
                                               "19"="A",
                                               "20"="J",
                                               "21"="H",
                                               "22"="H",
                                               "23"="E")
#gene projections
plot_cells(cds_subset, genes=c("zpg","aub","fzo","p-cup"), show_trajectory_graph = FALSE)

#generate rows for heatmap for single-nucleus data
genes <- c("r−l","shu","zpg","Rcd−1r","aub","mod","RpL27","Fkbp39","RpL27A","eEF1alpha1","Unr","RpS20","RpS20","CG32971","Nph","Pep","CG42857","RpS27A","RpL8","RpS25","RpL18A","RpS7","RpS8","RpL13A","RpS15","RpL6","smt3","vig2","RpL7","RpS11","RpS4","RpL19","RpL28","CG4415","Hel25E","sta","RpLP0","RpL10","RpS5b","RpLP1","RpL3","SC35","ATPsynC","tsr","ATPsynepsilonL","alphaTub84B","CG12493","Hsp26","Prosalpha6","Fkbp12","Prosbeta4","RpL37b","Ssb−c31a","Rbp4","Nlp","Prosbeta5","blanks","Prosbeta2","CycB","CG32259","Mtl","yps","MtnA","CG4021","CG42830","CG15876","fzo","CG5043","CycB","CG43800","CG31406","CG31093","CG15109","CG14113","CG42703","CG32137","CG17329","tHMG1","ATPsynCF6L","CG6652","CG14579","CG7768","CG18628","alphaTub84B","CG13747","CG42393","CG31523","CG17261","CG18662","CG12699","CG10934","CG42523","Pof","Jupiter","CG15657","COX4L","S−Lap2","Dic2","ProtA","CG31788","S−Lap3","CG14658","sowi","mil","CG31709","S−Lap8","CG31949","tHMG2","vrs","CG8840","Hsp60C","Mst33A","CG30393","CG4375","CG30039","CG5968","salto","knon","CG8838","CG17118","CG13898","CG13168","CG10859","CG5762","goddard","Mst77F","CG33293","betaTub85D","CG5280","S−Lap1","CG31948","CG5614","CG31294","CG14676","CG18418","CG15256","CG43183","CG13110","CG10734","ProtB","Ran−like","Mst84Da","CG12853","ocn","CG46385","CG31820","CG34168","CG46059","Mst84Dd","CG33340","CG31639","tomboy20","Mst98Ca","CG9129","CG9920","CG2127","CG14995","Pen","CG31226","CG3124","loopin−1","CG4691","CG10252","CG5089","CG17376","CG31988","CG17377","CG31538","CG17470","robl62A","CG15219","CG43935","CG31740","CG31468","CG12860","CG12861","dj","Hsp60B","S−Lap4","CG8701","janB","CG31407","eIF4E4","CG42688","CG30430","Mst84Db","Mst84Dc","Mst87F","CG9016","MCU","CG17666","CG18266","ste24c","CG8654","CG15025","CG17564","SdhAL","CG16984","CG1999","CG14605","Adgf−B","Rab3−GAP","CG31784","CG14689","m−cup","glob3","CG10749","CG14739","CG15638","CG8136","CG32161","CG13337","schuy","CG31730","CG16888","CG5906","CG11298","PCB","CG14505","CG30416","CG32228","CG32106","soti","loqs","CG11226","CG12448","CG30325","CG42570","CG12209","CG15548","boly","CG32396","CG11286","CG34167","CG43691","Pglym87","Tob","CG14864","CG42758","CG31817","CG30484","CG34216","CG13127","CG6652","c−cup","CG31735","sunz","CG31128","CG14106","wa−cup","d−cup","CG13538","CG15767","PpD6","p−cup","CG30178","Parp16","orb","CG30398","Glut3","CG30110","CG30334","CG8851","CG43796","CG14546","w−cup","CG6333","CG43339","CG31404","CG32681","Pif2","CG18368","CG10512","f−cup","CG15059","CG13330","CG30099","whip","CG12126","CG31815","CG42266","Efhc1.2","CG10969","CG33017","CG4073","r−cup","hale","CG31797","CG11629","CG43449","CG14294")
genes <- chartr("−","-",genes)

pt.matrix_nuc <- exprs(FCA_cds_subset)[match(genes,(rowData(FCA_cds_subset)$gene_short_name)),order(pseudotime(FCA_cds_subset))]
pt.matrix_nuc <- t(apply(pt.matrix_nuc,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_nuc <- t(apply(pt.matrix_nuc,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix_nuc) <- genes;

#generate annotation bar: nucleus
nuc_ann <- data.frame(FCA_cds_subset$cluster_id, pseudotime(FCA_cds_subset))
index <- order(nuc_ann$pseudotime.FCA_cds_subset.)
nuc_ann <- nuc_ann[index,]
nuc_ann <- data.frame(nuc_ann$FCA_cds_subset.cluster_id)
colours <- list('nuc_ann.FCA_cds_subset.cluster_id' = c('8'='#F8766D', 
                                                        '7'='#DB8E01', 
                                                        '6'='#AFA200',
                                                        '10'='#64B201',
                                                        '11'='#64B201',
                                                        '13'='#02BD5D',
                                                        '12'='#02C1A7',
                                                        '4' ='#02C1A7',
                                                        '5'='#02C1A7',
                                                        '9'='#02C1A7',
                                                        '2'='#02BADE',
                                                        '1'='#00A6FF',
                                                        '14'='#EF67EB',
                                                        '3'='#FF62B6'))
colAnn_nuc <- HeatmapAnnotation(df = nuc_ann, which = 'col', col=colours, annotation_width = unit(c(1, 4), 'cm'), gap = unit(1, 'mm'))
htkm_adult <- Heatmap(
  pt.matrix_nuc,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 5,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = colAnn_nuc)

htkm_adult

#generate rows for heatmap for single-cells data
pt.matrix_cell <- exprs(cds_subset)[match(genes,(rowData(cds_subset)$gene_short_name)),order(pseudotime(FCA_cds_subset))]
pt.matrix_cell <- t(apply(pt.matrix_nuc,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix_cell <- t(apply(pt.matrix_nuc,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix_cell) <- genes;

#annotation bar
ann <- data.frame(cds_subset$cluster_id, pseudotime(cds_subset))
index <- order(ann$pseudotime.cds_subset.)
ann <- ann[index,]
ann <- data.frame(ann$cds_subset.cluster_id)
colours <- list('ann.cds_subset.cluster_id' = c('A' = '#F8766D', 
                                                'B' = '#DB8E01', 
                                                'C'='#AFA200',
                                                'D'='#64B201',
                                                'E'='#02BD5D',
                                                'F'='#02C1A7',
                                                'G'='#02BADE',
                                                'F'='#02C1A7',
                                                'G'='#02BADE',
                                                'H'='#00A6FF',
                                                'I'='#B385FF',
                                                'J'='#EF67EB',
                                                'K'='#FF62B6'))
colAnn <- HeatmapAnnotation(df = ann, which = 'col',col=colours,annotation_width = unit(c(1, 4), 'cm'), gap = unit(1, 'mm'))

#heatmap for cell data
htkm_cell <- ComplexHeatmap::Heatmap(
  pt.matrix_cell,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 1,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE, 
  top_annotation               = colAnn)

htkm_nuc
htkm_cell

#ggplot of warped pseudotime (from SP)
library(patchwork, lib.loc = "/nfs/apps/lib/R/3.6-focal/site-library.2021q2")

df <- readRDS('/lab/solexa_yamashita/people/Amelie/FCA_data/warped_data.rds')

type1_genes = c('S-Lap2', 'Dic2', 'ProtA', 'CG31788','S-Lap3','CG14658', 'sowi')

type2_genes = c('CG30430', 'Mst84Db', 'Mst84Dc', 'Mst87F', 'CG9016', 'MCU', 'CG17666')

get_plot <- function(df, gene_list) {
  (
    df
    %>% dplyr::filter(gene %in% gene_list)
    %>% gather('pseudotime', 'expression', -gene, -DataSet)
    %>% mutate(pseudotime = as.numeric(pseudotime))
    %>% mutate(label = factor(gene, levels=gene_list, ordered=TRUE))
    %>% mutate(DataSet = factor(DataSet, levels=c("AdultCell", "AdultNuclei"), ordered=TRUE))
    %>% ggplot(aes(pseudotime, expression, color=DataSet))
    + geom_line(size=1)
    + facet_wrap(~label, ncol=7)
    + labs(x="warped pseudotime", y="z-score")
    + scale_color_manual(values = c("AdultCell"="#0000CC", "AdultNuclei"="#CC0000"))
    + theme_minimal()
    + theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(colour = "white", fill=NA)
    )
  )
}

p1 <- get_plot(df, type1_genes)
p2 <- get_plot(df, type2_genes)

p <- p1 + p2 + plot_layout(ncol=1, nrow = 2) & theme(legend.position="right")
p

ggsave("/lab/solexa_yamashita/people/Amelie/FCA_data/selected_genes_warp_pattern_two.pdf", p, width=8.4, height=3.5)



