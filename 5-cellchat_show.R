rm(list = ls())
library(Seurat)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)

load("1-cellchat.rda")
df.net <- subsetCommunication(cellchat)

#### 1.Systems analysis of cell-cell communication network
cellchat <- netAnalysis_computeCentrality(cellchat)

netAnalysis_signalingRole_scatter(cellchat,group = c("TAM","DC"))

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

table(cellchat@netP$pathways)
table(cellchat@idents)


pathways.show <- c("ADGRE5","ALCAM","ANNEXIN","APP","CCL","CD137","CD45","CD6","CD86","CD99",
                   "CXCL","GALECTIN","ICAM","IFN-II","IL1","IL16","ITGB2","JAM","LCK",
                   "MIF","NEGR","NOTCH","PECAM1","SELPLG","SEMA4","SN","TGFb","THBS",
                   "TNF","VCAM","VEGF","VISFATIN")#"MHC-I","MHC-II",

pathways.show <- cellchat@netP$pathways
pathways.show <- pathways.show[pathways.show!="MHC-II"]
netVisual_bubble(cellchat, sources.use = c("PD1+TRM","MAIT","IFNG+Tc")
                 , targets.use = c("TAM"),signaling = pathways.show)

pathways.show <- cellchat@netP$pathways
pathways.show <- pathways.show[pathways.show!="MHC-I"]
netVisual_bubble(cellchat, targets.use = c("PD1+TRM","MAIT","IFNG+Tc")
                 , sources.use = c("TAM"),signaling = pathways.show)