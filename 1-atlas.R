rm(list = ls())
library(dplyr)
library(Seurat)
library(data.table)

hcc.list <- c(cohort1,cohort2,cohort3)

for (i in 1:length(hcc.list)) {
  hcc.list[[i]] <- subset(hcc.list[[i]], subset = nFeature_RNA > 500 & percent.mt < 25)
}

hcc.list

hcc.anchors <- FindIntegrationAnchors(object.list = hcc.list, dims = 1:30)
hcc.integrated <- IntegrateData(anchorset = hcc.anchors, dims = 1:30)

DefaultAssay(hcc.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
hcc.integrated <- ScaleData(hcc.integrated, verbose = FALSE)
hcc.integrated <- RunPCA(hcc.integrated, npcs = 30, verbose = FALSE)

DimPlot(hcc.integrated, reduction = "pca",group.by = "COHORT")

hcc.integrated <- FindNeighbors(hcc.integrated, dims = 1:10)
hcc.integrated <- FindClusters(hcc.integrated, resolution = 0.5)
table(hcc.integrated$seurat_clusters)

hcc.integrated <- RunUMAP(hcc.integrated, reduction = "pca", dims = 1:30)

DimPlot(hcc.integrated, reduction = "umap",label = T)
DimPlot(hcc.integrated, reduction = "umap",group.by = "COHORT")
DimPlot(hcc.integrated, reduction = "umap",group.by = "SAMPLE")
DimPlot(hcc.integrated, reduction = "umap",group.by = "TISSUE")
DimPlot(hcc.integrated, reduction = "umap",group.by = "Type")
DimPlot(hcc.integrated, reduction = "umap",group.by = "SAMPLE")

total <- hcc.integrated

library(SingleR)

load("../package/SingleR/hpca.rda")
refdata <- hpca.se

testdata <- GetAssayData(total, slot="data")
clusters <- total@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

total@meta.data$celltype_hpca = "NA"
for(i in 1:nrow(celltype)){
  total@meta.data[which(total@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_hpca'] <- celltype$celltype[i]}

DimPlot(total, group.by="celltype_hpca", repel=T, label=T, reduction='umap')

save(total,file = "1-atlas.rda")