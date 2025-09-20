rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)

##############
load("CD8T.rda")
liver <- CD8T
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
DimPlot(liver, reduction = "umap",group.by = "SAMPLE")

total.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(total.markers,file = "1-totalmarkers.rda")

DotPlot(liver, features = c("CD69","ITGAE","CXCR6","ITGA1","GZMB","IFNG","IL17A",
                            "CTLA4","PDCD1","TIGIT","PCNA"), assay = "RNA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_gradient(low = "white", high = "#E88B20")

liver$MM3 <- "NA"
liver$MM3[liver$seurat_clusters%in%"0"] <- "GZMH+Tc"
liver$MM3[liver$seurat_clusters%in%"1"] <- "TEM"
liver$MM3[liver$seurat_clusters%in%"2"] <- "IFNG+Tc"
liver$MM3[liver$seurat_clusters%in%"3"] <- "RGS1+TRM"
liver$MM3[liver$seurat_clusters%in%"4"] <- "CCL4+Tc"
liver$MM3[liver$seurat_clusters%in%"5"] <- "MAIT"
liver$MM3[liver$seurat_clusters%in%"6"] <- "PD1+TRM"
liver$MM3[liver$seurat_clusters%in%"7"] <- "TEX"
liver$MM3[liver$seurat_clusters%in%"8"] <- "XCL2+Tc"
DimPlot(liver, reduction = "umap",group.by = "MM3",label = T)

save(liver,file = "1-CD8T.rda")
