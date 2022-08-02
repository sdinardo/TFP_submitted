# Purpose:  Extract Cyst lineage, re-run UMAP to enable 3D plotting 
# 2021-20-07 The R script, cystLineage_in3D_2021_10.R, was written by SDFlies with help from @Amelie (code for re-running UMAP)

# The cystLineage extracted below uses @Jennifer Viveiros's suggestion of 22 leiden 6 clusters restricted to the 'true' cystLineage;  using those clusters, the cystLineage contains both 'feet' and now excludes the mermaid tail (it excludes leiden 6 clusters 20, 42, 89).

#2022_07: NOTE --> This object can now be created by reading in an updated testis FCA SEURAT object: FCAloomToSeurat2TFP_Annotations.rds; 
#  You may need to change fileNames and or filePaths in the code below.

#  This is used for Figure 6C.

# To run the script, you need the testis FCA SEURAT object in your environment. Either readRDS 'FCA_bioHub_ScrapeB.rds', or the newer 'FCA_79hubcells.rds' (made by @Gabriela Vida), which separates old Cluster 90 (of 120 cells) into:
# 1. True hub of 79 cells, still called Cluster 90;
# 2. Plus 41 extraneous cells, now labeled Cluster 111.
# Ask Gabriela Vida (GV) for access to a Box folder with these .rds files.
# https://upenn.app.box.com/folder/145907944818?s=zc0zxqc2p81twzw034kp34q81vu490fh




# If you just wish to interrogate the cystLineage object, the script below produces an .rds file you can also ask for (instead of producing yourself):
# cystLineage_for3D_2021_10.  

# A 2021_10_29  update to the script fixed a problem; now the 3D UMAP parallels the original 3D map in the full FCA project. 


# Load packages   ####
library(Seurat)
library(tidyverse)
library(plotly)
library(ggplot2)

FCA_79hubcells <- readRDS(file = "~/Desktop/FCA R/FCA_79hubcells.rds")

# Subset; isolate the Cyst lineage   #### 
cystLineage2021_10 <- subset(FCA_79hubcells, idents = c("62", "36", "58", "77", "65", "47","95","88","57","80","84","72","74","98","69","67","60","94","104","81","79","68")) # 22 clusters, restricted to cystLineage  containing both feet, and excluding mermaid tail (20, 42, 89). 

cystLineage2021_10[["leiden6.ident"]] <- Idents(object = cystLineage2021_10) # Store cluster identities to reuse in the subsetted SEURAT object.

# Run a standard analysis pipeline   ####
# The subsetted object inherits both "variable Featerus and scaled dfata fromthe original object.
# Thus we should skip these two steps:
# cystLineage2021_10 <- FindVariableFeatures(cystLineage2021_10, selection.method = "vst", nfeatures = 2000)
# cystLineage2021_10 <- ScaleData(cystLineage2021_10, features = all.genes)

cystLineage2021_10 <- RunPCA(cystLineage2021_10, features = VariableFeatures(object = cystLineage2021_10))
cystLineage2021_10 <- FindNeighbors(cystLineage2021_10, dims = 1:10)
cystLineage2021_10 <- FindClusters(cystLineage2021_10, resolution = 0.5)

# Runa a new UMAP, passing in an argument to generate 3 dimensions instead of only 2: 
cystLineage2021_10 <- RunUMAP(cystLineage2021_10, dims = 1:10, n.components = 3L) # n.components = 3L' will produce all 3 dimensions for true 3D plotting.

Idents(object = cystLineage2021_10) <- "leiden6.ident" # re-establish our original cluster identities 

colorDimPlot <- DimPlot(cystLineage2021_10, label = TRUE) # Check the plot
DimPlot(cystLineage2021_10, label = TRUE) 
ggsave("cystlineage2D_3Dcode.jpg", plot = last_plot(), device = jpeg, scale = 1, dpi = 600)
# Save the new object   ####
# Uncomment this next line to save the object to file
# Change the path below
saveRDS(cystLineage2021_10, file = "~/Deskop/FCA R/cyst lineage 3D mapping 12.14/cystLineage_for3D_2021_10.rds")

# Construct a dataframe to use in plotting   ####
# We will pull various bits of info from our subsetted SEURAT object in order to produce the requisite plots
# 1. UMAP coordinates
# 2. Expression data for a couple markers (as examples)
# 3. The cluster numbers


# 1.  Extract umap coordinates in order to plot in 3D

# Visualize what the headings are called so that you can extract them to form the dataframe
head(Embeddings(object = cystLineage2021_10, reduction = "umap")) # prints cellIDs, and a column for each UMAP dimension's coordinates.

# Extract the UMAP 'embeddings' from the Seurat Object
umap_1 <- cystLineage2021_10[["umap"]]@cell.embeddings[,1]  # Asking for 'all rows' in column 1
umap_2 <- cystLineage2021_10[["umap"]]@cell.embeddings[,2]
umap_3 <- cystLineage2021_10[["umap"]]@cell.embeddings[,3]

# Prepare a dataframe for cell plotting

# 2.  Create a vector of the example genes of interest (goi):
goi <- c("tj", "so")
# Use 'FetchData' tp pull various columns:
plotting.data <- FetchData(object = cystLineage2021_10, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi, "leiden6.ident"), slot = 'data') # 'slot' indicates which type of expression data to pull. "data" = normalized; "scale.data" = norm+scaled.

head(plotting.data) # Inspect the dataframe, to understand what columns you will be plotting.

# Plot a 3D UMAP   ####
# We will use Ployly, and it would be convenient to keep 'colors' the same were used in the DimPlot: https://github.com/satijalab/seurat/issues/257

# Load the "scales" package
require(scales)

# Create a vector with levels of object@ident
identities <- levels(cystLineage2021_10@meta.data$leiden6.ident)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

fig <- plot_ly(data = plotting.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~leiden6.ident,
               colors = my_color_palette, # we created this to match DimPot colors
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 2, width=2), # controls the size of points
               text=~leiden6.ident, # Specifies what will show upon 'hover'
               hoverinfo="text") # Hovering your mouse will show cluster names

fig # Prints the figure; to save --> Export > Save as Web Page.


# Next, plot gene expression in 3D   ####

# First check how the usual 2D plot looks: 
FeaturePlot(cystLineage2021_10, features = goi, pt.size = 0.8, ncol = 2)  # usual 2D plot of gene expression


# Using plot_ly, check a 3D plot of one gene's expression:
figCyst_tj <- plot_ly(data = plotting.data,
                      x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                      color = ~tj,
                      opacity = 0.9,
                      colors = c('lightgrey', 'darkblue'),
                      type = "scatter3d", 
                      mode = "markers",
                      marker = list(size = 3), 
                      text=~leiden6.ident,
                      hoverinfo="text"
) #  %>%layout(title= "tj & so")


# Using plot_ly, check out a 3D plot of a second gene's expression:
figCyst_so <- plot_ly(data = plotting.data,
                      x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                      color = ~so,
                      opacity = 0.3,
                      colors = c('lightgrey', 'red'),
                      type = "scatter3d", 
                      mode = "markers",
                      marker = list(size = 2), 
                      text=~leiden6.ident,
                      hoverinfo="text"
) #  %>%layout(title= "tj & so")


# Now superimpose the two plots
figTjSo <- subplot(figCyst_tj, figCyst_so) %>% 
  layout(title = 'Tj & So in Cyst lineage')

figTjSo  # Prints the figure; to save --> Export > Save as Web Page.

# Note: there are likely better ways to map multiple genes on the 3D plot.
© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
