rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)


##############
load("Myeloid_cell.rda")
liver <- Myeloid_cell
table(liver$SAMPLE)

liver <- RunPCA(liver, features = VariableFeatures(object = liver))
DimPlot(liver, reduction = "pca")
DimPlot(liver, reduction = "pca",group.by = "TISSUE")


liver <- FindNeighbors(liver, dims = 1:10)
liver <- FindClusters(liver, resolution = 0.5)

table(liver$seurat_clusters)

liver <- RunUMAP(liver, dims = 1:10)

DimPlot(liver, reduction = "umap",label = T)
DimPlot(liver, reduction = "umap",group.by = "TISSUE")
DimPlot(liver, reduction = "umap",group.by = "COHORT")
DimPlot(liver, reduction = "umap",group.by = "TREATMENT2")
DimPlot(liver, reduction = "umap",group.by = "MM3",label = T)


liver$MM3 <- "NA"
liver$MM3[liver$seurat_clusters%in%"0"] <- "TAM1"
liver$MM3[liver$seurat_clusters%in%"1"] <- "TAM2"
liver$MM3[liver$seurat_clusters%in%"2"] <- "TAM3"
liver$MM3[liver$seurat_clusters%in%"3"] <- "TAM4"
liver$MM3[liver$seurat_clusters%in%"4"] <- "DC1"
liver$MM3[liver$seurat_clusters%in%"5"] <- "TAM5"
liver$MM3[liver$seurat_clusters%in%"6"] <- "TAM6"
liver$MM3[liver$seurat_clusters%in%"7"] <- "TAM7"
liver$MM3[liver$seurat_clusters%in%"8"] <- "DC2"
liver$MM3[liver$seurat_clusters%in%"9"] <- "TAM8"
liver$MM3[liver$seurat_clusters%in%"10"] <- "pDC"
liver$MM3[liver$seurat_clusters%in%"11"] <- "DC3"
DimPlot(liver, group.by="MM3", repel=T, reduction='umap',label = T)

save(liver,file = "1-atlas.rda")
