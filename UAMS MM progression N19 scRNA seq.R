library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(pheatmap)
library(pals)
library(ggplot2)
library(devtools)
library(fields)
library(MAST)
library(KernSmooth)


#################################################
###### Load and QC data
#################################################
options(future.globals.maxSize = 100000 * 1024^2)
# HD_1
import <- Read10X_h5("HD_1.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
HD_1 <- import
rm(import)

# HD_2
import <- Read10X_h5("HD_2.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
HD_6 <- import
rm(import)

# HD_3
import <- Read10X_h5("HD_3.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
HD_3 <- import
rm(import)

# HD_4
import <- Read10X_h5("HD_4.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
HD_4 <- import
rm(import)

# HD_5
import <- Read10X_h5("HD_5.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
HD_5 <- import
rm(import)

# MGUS_1
import <- Read10X_h5("MGUS_1.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MGUS_1 <- import
rm(import)

# MGUS_2
import <- Read10X_h5("MGUS_2.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MGUS_2 <- import
rm(import)

# MGUS_3
import <- Read10X_h5("MGUS_3.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MGUS_9 <- import
rm(import)

# MGUS_4
import <- Read10X_h5("MGUS_4.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MGUS_10 <- import
rm(import)

# MGUS_5
import <- Read10X_h5("MGUS_5.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MGUS_11 <- import
rm(import)

# MGUS_6
import <- Read10X_h5("MGUS_6.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MGUS_13 <- import
rm(import)

# SMM_1
import <- Read10X_h5("SMM_1.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
SMM_1 <- import
rm(import)

# SMM_2
import <- Read10X_h5("SMM_2.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
SMM_2 <- import
rm(import)

# SMM_3
import <- Read10X_h5("SMM_3.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
SMM_3 <- import
rm(import)

# SMM_4
import <- Read10X_h5("SMM_4.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
SMM_5 <- import
rm(import)

# MM_1
import <- Read10X_h5("MM_1.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MM_1 <- import
rm(import)

# MM_2
import <- Read10X_h5("MM_2.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MM_3 <- import
rm(import)

# MM_3
import <- Read10X_h5("MM_3.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MM_3 <- import
rm(import)

# MM_4
import <- Read10X_h5("MM_4.h5")
import <- CreateSeuratObject(counts = import, assay = 'RNA', project = 'MGUS', min.features = 200)
import[["percent.mt"]] <- PercentageFeatureSet(import, pattern = "^mt-")
import <- subset(import, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)
MM_4 <- import
rm(import)

####################
#### integrate data
####################

# list of Seurat object
cart.list <- list(MGUS_1, MGUS_2, MGUS_3, MGUS_4, MGUS_5, MGUS_6,
                  SMM_1, SMM_2, SMM_3, SMM_4,
                  MM_1, MM_2, MM_3, MM_4,
                  HD_1, HD_2, HD_3, HD_4, HD_5)

# normalize and identify variable features for each dataset independently
cart.list <- lapply(X = cart.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

## fasta integration using rpca
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = cart.list)
cart.list <- lapply(X = cart.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#################################################
###### Perform integration
#################################################

cart.anchors <- FindIntegrationAnchors(object.list = cart.list, anchor.features = features, reduction = "rpca")

cart.combined <- IntegrateData(anchorset = cart.anchors)
DefaultAssay(cart.combined) <- "integrated"

cart.combined <- ScaleData(cart.combined, verbose = FALSE)
cart.combined <- RunPCA(cart.combined, npcs = 50, verbose = FALSE)
cart.combined <- FindNeighbors(cart.combined, reduction = "pca", dims = 1:30)
cart.combined <- FindClusters(cart.combined, resolution = 0.2)
cart.combined <- RunUMAP(cart.combined, reduction = "pca", dims = 1:20, min.dist = 0.5, n.neighbors = 30)

Idents(cart.combined) <- 'state'
Idents(cart.combined) <- factor(Idents(cart.combined), levels = c('HD','MGUS', 'SMM', 'MM'))
cart.combined$state <- factor(x = cart.combined$state, levels = c('HD','MGUS', 'SMM', 'MM'))
levels(cart.combined)
Idents(cart.combined) <- 'seurat_clusters'

DimPlot(cart.combined, reduction = "umap", label = T, label.size = 4,  repel = F)

DimPlot(cart.combined, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)


####################
#### Remove Erythrocyte
####################

New_WOErythrocyte <- subset(cart.combined, idents = c('0','1','2','4','5','7','10','12','13','14','15','16','17','18','20','21','22'))
New_WOErythrocyte <- RunPCA(New_WOErythrocyte, npcs = 50, verbose = FALSE)
New_WOErythrocyte <- FindNeighbors(New_WOErythrocyte, reduction = "pca", dims = 1:50)
New_WOErythrocyte <- FindClusters(New_WOErythrocyte, graph.name = 'integrated_snn', resolution = 0.2) # resolution: number of clusters
New_WOErythrocyte <- RunUMAP(New_WOErythrocyte, reduction = "pca", dims = 1:50, min.dist = 0.9, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(New_WOErythrocyte, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(New_WOErythrocyte, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)

####################
#### Final_Clus
####################

Final_Clus <- subset(New_WOErythrocyte, idents = c('0','1','2','3','4','5','6','8','9','10','11','12','13','14','15','16','18'))
Final_Clus <- RunPCA(Final_Clus, npcs = 50, verbose = FALSE)
Final_Clus <- FindNeighbors(Final_Clus, reduction = "pca", dims = 1:50)
Final_Clus <- FindClusters(Final_Clus, graph.name = 'integrated_snn', resolution = 0.2) # resolution: number of clusters
Final_Clus <- RunUMAP(Final_Clus, reduction = "pca", dims = 1:50, min.dist = 0.9, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(Final_Clus, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(Final_Clus, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)


####################
#### combine clusters
####################

# change all ids to the same cluster id : ('1','4','25) -> ('1','1','1')

new.cluster.ids <- c('0','1','1','2','3','4','1','5','6','7','8','9','1','3','10','11','0')

names(new.cluster.ids) <- levels(Final_Clus)
Final_Clus <- RenameIdents(Final_Clus, new.cluster.ids)

Final_Clus$seurat_clusters <- Idents(Final_Clus)

DimPlot(Final_Clus, reduction = "umap", label = T, label.size = 4,  repel = F)

DimPlot(Final_Clus, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)

DimPlot(Final_Clus, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 4)

rm(cart.list, MGUS_1, MGUS_2, MGUS_9, MGUS_10, MGUS_11, MGUS_13,
   SMM_1, SMM_2, SMM_3, SMM_5,
   MM_1, MM_3, MM_4, MM_6,
   HD_1, HD_3, HD_4, HD_5, HD_6)

rm(cart.anchors, cart.list)
rm(New_WOErythrocyte)
rm(cart.combined)

##DotPlot
Final_Clus$sample_id <- factor(Final_Clus$sample_id, levels = c('MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 
                                                              'MGUS_5', 'SMM_1',
                                                              'MGUS_6', 'MM_1',
                                                              'SMM_2', 'MM_2',
                                                              'SMM_3', 'MM_3',
                                                              'SMM_4', 'MM_4',
                                                              'HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5'))
Idents(Final_Clus) <- 'seurat_clusters'
Idents(Final_Clus)


DotPlot(Final_Clus, assay = 'RNA', features = c('ELK1', 'SF3B1', 'CEP55', 'NEK2', 'TOP2A', 'CHML', 'RRM2', 'PBK', 'TRIP13', 'KIF4A', 'CDC20'), 
        cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )


DefaultAssay(Final_Clus) <- 'RNA'
FeaturePlot(Final_Clus, features = c('MX1'), split.by = 'state', min.cutoff = 'q25') & theme(text = element_text(size=10),
                                                                          axis.title = element_text(size=12),
                                                                          legend.text=element_text(size=12),
                                                                          legend.title=element_text(size=12)
)


##cluste all.markers

Final_Clus <- JoinLayers(Final_Clus, assay = 'RNA')
all.markers <- FindAllMarkers(object = Final_Clus, test.use = 'MAST', min.pct = 0.2, 
                              min.cells.feature = 10,
                              min.cells.group = 10,
                              only.pos = TRUE)


top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, file = 'Final_Clus All DEG.csv', row.names = FALSE)


#################################################
###### Subcluster analysis
#################################################

##
# MM cells
##

MM_cluster <- subset(Final_Clus, idents = c('2'))
MM_cluster <- RunPCA(MM_cluster, npcs = 50, verbose = FALSE)
MM_cluster <- FindNeighbors(MM_cluster, reduction = "pca", dims = 1:50)
MM_cluster <- FindClusters(MM_cluster, graph.name = 'integrated_snn', resolution = 2.5) # resolution: number of clusters
MM_cluster <- RunUMAP(MM_cluster, reduction = "pca", dims = 1:50, min.dist = 0.5, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(MM_cluster, reduction = "umap",label = T, label.size = 4, repel = F)

DimPlot(MM_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 5)

DimPlot(MM_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)


####################
#### combine clusters 2nd
####################

# change all ids to the same cluster id : ('1','4','25) -> ('1','1','1')

MM_new.cluster.ids <- c('0','1','2','2','3','4','5','6','4','7',
                        '8','9','10','4','11','12','12','13','12','14','12',
                        '15','16','8','5')

names(MM_new.cluster.ids) <- levels(MM_cluster)
MM_cluster <- RenameIdents(MM_cluster, MM_new.cluster.ids)

MM_cluster$seurat_clusters <- Idents(MM_cluster)

DimPlot(MM_cluster, reduction = "umap", label = T, label.size = 4,  repel = F)

DimPlot(MM_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 4)

###
# Sub cluste all.markers
###
Idents(MM_cluster) <- 'seurat_clusters'
Idents(MM_cluster)
MM_cluster <- JoinLayers(MM_cluster, assay = 'RNA')
all.markers <- FindAllMarkers(object = MM_cluster, test.use = 'MAST', min.pct = 0.2, 
                              min.cells.feature = 10,
                              min.cells.group = 10,
                              only.pos = TRUE)


top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, file = 'MM_cluster17_All.csv', row.names = FALSE)

##################
### subset to identify different group
##################

MM_cluster_subset <- subset(MM_cluster, idents = c('HD','MGUS','SMM'))

##########
# find DEG between groups
##########

Idents(MM_cluster) <- 'state'
Idents(MM_cluster)

## markers
MM_cluster <- JoinLayers(MM_cluster, assay = 'RNA')
DEG_MGUSvsHD <- FindMarkers(MM_cluster, assay = 'RNA', slot = 'data',
                            ident.1 = 'HD',
                            test.use = 'MAST',
                            min.pct = 0.2, only.pos = FALSE, verbose = TRUE)

DEG_MGUSvsHD$gene <- rownames(DEG_MGUSvsHD)

write.csv(DEG_MGUSvsHD, file = 'MM_cluster_DEG_MGUSvsHD.csv', row.names = FALSE) ## write CSV

## Violin plot
DefaultAssay(MM_cluster) <- 'RNA'
VlnPlot(MM_cluster, features = c("PCDH9"), pt.size = 0)

##Sub cluster Counting
rep_cell_counts <- group_by(MM_cluster@meta.data, seurat_clusters, sample_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "MM_cluster 17 Cell counts.csv",sep = ",",row.names = FALSE)

## Feature plot "butterfly plot"
DefaultAssay(MM_cluster) <- 'RNA'
FeaturePlot(MM_cluster, features = c('CD37')) & theme(text = element_text(size=10),
                                                     axis.title = element_text(size=12),
                                                     legend.text=element_text(size=12),
                                                     legend.title=element_text(size=12)
)

###
# Sub cluster top20.markers
###
Idents(MM_cluster) <- 'seurat_clusters'
Idents(MM_cluster)

all.markers <- FindAllMarkers(object = MM_cluster, test.use = 'MAST', min.pct = 0.2, only.pos = TRUE)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(top20.markers, file = 'MM_cluster_TOP100.csv', row.names = FALSE)

 
##Dot Plot by sample_id
MM_cluster$sample_id <- factor(MM_cluster$sample_id, levels = c('MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 
                                                                'MGUS_5', 'SMM_1',
                                                                'MGUS_6', 'MM_1',
                                                                'SMM_2', 'MM_2',
                                                                'SMM_3', 'MM_3',
                                                                'SMM_4', 'MM_4',
                                                                'HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5'))
Idents(MM_cluster) <- 'sample_id'


DotPlot(MM_cluster, assay = 'RNA', features = c('EDNRB', 'CCND2', 'CCND1', 'DKK1', 'FRZB', 'CNTN5', 'NCAM1', 'PCDH9', 'CADM1', 'IFI27', 'IFIT1',
                                                'IFI6', 'IFITM1', 'ISG15', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DRB1', 'HSPA1B', 'HSPA1A', 'IGKV4-1',
                                                'IGKV3-11', 'IGLC3', 'IGLC2', 'IGLC1', 'IGLV6-57', 'IGLV3-1', 'IGHV1-69D', 'IGHA2', 'IGHA1',
                                                'IGHD', 'IGHG4', 'IGHG3', 'IGHG2', 'IGHG1', 'IGHGP', 'IGHM'), 
        cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

##Dot Plot by cluster

Idents(MM_cluster) <- 'seurat_clusters'
Idents(MM_cluster)

DotPlot(MM_cluster, assay = 'RNA', features = c('IGHD', 'PCDH9', 'IGHV1-69D', 'RPS4Y1', 'VPREB3', 'EDNRB', 'LAMP5', 'LPAR6', 'DKK1', 'IGKV3-11',
                                                'NUCB2', 'SSR4', 'CD63', 'FKBP11', 'DERL3', 'EIF1', 'SSR3', 'PSAT1', 'SERP1', 'SEC61B',
                                                'IGHV4-59', 'IGHG3', 'IGHG4', 'IGHG1', 'DPEP1', 'SNRPG', 'MRPL52', 'COX14', 'EDEM2', 'GLRX5',
                                                'HEXIM1', 'Z93241.1', 'MT-ND6', 'ZFAS1', 'INTS6', 'IGKV1D-16', 'CCND1', 'IGKV1-16', 'MS4A1', 'SERPINE2',
                                                'NEB', 'IGHV4-34', 'CST6', 'IGKV4-1', 'IGHG2', 'ZNF91', 'RAB30-DT', 'DDX17', 'SNHG14', 'ADA2',
                                                'RRBP1', 'XIST', 'MT-ATP8', 'MEF2C', 'MT-CO2', 'CD74', 'LINC01781', 'CD37', 'HLA-DRA', 'HLA-E',
                                                'HNRNPC', 'DDX24', 'TPD52', 'UBXN4', 'TOP1', 'IGHGP', 'CCND2', 'LINC01229', 'PRDX1', 'UAP1',
                                                'CCR10', 'CSDE1', 'ME2', 'BIRC3', 'GSTP1', 'LINC02432', 'MAP3K14', 'HSP90B1', 'FTH1', 'PCGF5',
                                                'WFDC2', 'RBP1', 'SPINK2'), 
        cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text()  # Rotate x-axis labels by 45 degrees
  )


Idents(MM_cluster) <- 'state'
DefaultAssay(MM_cluster) <- 'RNA'
VlnPlot(MM_cluster, features = c("CC"), split.by = 'state')

##
# ImG cells
##

ImG_cluster <- subset(Final_Clus, idents = c('5'))
ImG_cluster <- RunPCA(ImG_cluster, npcs = 50, verbose = FALSE)
ImG_cluster <- FindNeighbors(ImG_cluster, reduction = "pca", dims = 1:50)
ImG_cluster <- FindClusters(ImG_cluster, graph.name = 'integrated_snn', resolution = 0.2) # resolution: number of clusters
ImG_cluster <- RunUMAP(ImG_cluster, reduction = "pca", dims = 1:50, min.dist = 0.5, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(ImG_cluster, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(ImG_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)
DimPlot(ImG_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 5)

##########
# find DEG
##########

Idents(ImG_cluster) <- 'state'
Idents(ImG_cluster)

## markers
ImG_cluster <- JoinLayers(ImG_cluster, assay = 'RNA')
DEG_MGUSvsHD <- FindMarkers(ImG_cluster, assay = 'RNA', slot = 'data',
                            ident.1 = 'HD',
                            test.use = 'MAST',
                            min.pct = 0.2, only.pos = FALSE, verbose = TRUE)

DEG_MGUSvsHD$gene <- rownames(DEG_MGUSvsHD)

write.csv(DEG_MGUSvsHD, file = 'DEG_MGUSvsHD.csv', row.names = FALSE) ## write CSV

Idents(ImG_cluster) <- 'state'
DefaultAssay(ImG_cluster) <- 'RNA'
VlnPlot(ImG_cluster, features = c("ISG15"), split.by = 'state')


###
# Sub cluste top20.markers
###
Idents(ImG_cluster) <- 'seurat_clusters'
Idents(ImG_cluster)

all.markers <- FindAllMarkers(object = ImG_cluster, test.use = 'MAST', min.pct = 0.2, only.pos = TRUE)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
write.csv(top20.markers, file = 'ImG_cluster_TOP15.csv', row.names = FALSE)

##Sub cluster Counting
rep_cell_counts <- group_by(ImG_cluster@meta.data, seurat_clusters, sample_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "ImG_cluster Cell counts.csv",sep = ",",row.names = FALSE)


##DotPlot
ImG_cluster$sample_id <- factor(ImG_cluster$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                  'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                  'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                  'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(ImG_cluster) <- 'sample_id'


DotPlot(ImG_cluster, assay = 'RNA', features = c('ELK1', 'SF3B1', 'CEP55', 'NEK2', 'TOP2A', 'CHML', 'RRM2', 'PBK', 'TRIP13', 'KIF4A', 'CDC20'), 
        cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )


######
# T cells
######

T_cluster <- subset(Final_Clus, idents = c('1'))
T_cluster <- RunPCA(T_cluster, npcs = 50, verbose = FALSE)
T_cluster <- FindNeighbors(T_cluster, reduction = "pca", dims = 1:50)
T_cluster <- FindClusters(T_cluster, graph.name = 'integrated_snn', resolution = 0.6) # resolution: number of clusters
T_cluster <- RunUMAP(T_cluster, reduction = "pca", dims = 1:50, min.dist = 1.2, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(T_cluster, reduction = "umap",label = T, label.size = 4, repel = F) ##show UMAP plot
DimPlot(T_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)
DimPlot(T_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 5)



####################
#### combine clusters
####################

# change all ids to the same cluster id : ('1','4','25) -> ('1','1','1')

T_new.cluster.ids <- c('0','1','2','3','3','4','5','2','6','7','8','8',
                       '9','10','11','10','1','12','12','13','8','3','8','14')

names(T_new.cluster.ids) <- levels(T_cluster)
T_cluster <- RenameIdents(T_cluster, T_new.cluster.ids)

T_cluster$seurat_clusters <- Idents(T_cluster)

DimPlot(T_cluster, reduction = "umap", label = T, label.size = 4,  repel = F)

DimPlot(T_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 4)

DefaultAssay(T_cluster) <- 'RNA'
FeaturePlot(T_cluster, features = c('STAT1'),split.by = 'state') & theme(text = element_text(size=10),
                                                      axis.title = element_text(size=12),
                                                      legend.text=element_text(size=12),
                                                      legend.title=element_text(size=12)
)


T_cluster$sample_id <- factor(T_cluster$sample_id, levels = c('MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 
                                                                'MGUS_5', 'SMM_1',
                                                                'MGUS_6', 'MM_1',
                                                                'SMM_2', 'MM_2',
                                                                'SMM_3', 'MM_3',
                                                                'SMM_4', 'MM_4',
                                                                'HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5'))


###
# Sub cluster top20.markers
###
Idents(T_cluster) <- 'seurat_clusters'
Idents(T_cluster)

T_cluster <- JoinLayers(T_cluster, assay = 'RNA')

all.markers <- FindAllMarkers(object = T_cluster, test.use = 'MAST', min.pct = 0.2, only.pos = TRUE)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20.markers, file = 'Final_T_clusterr_TOP20.csv', row.names = FALSE)


##Sub cluster Counting
rep_cell_counts <- group_by(T_cluster@meta.data, seurat_clusters, sample_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "sub-cluster T Cell counts.csv",sep = ",",row.names = FALSE)


##Dot Plot
T_cluster$seurat_clusters <- factor(T_cluster$seurat_clusters, levels = c('3', '8', '10', '12', 
                                                                          '7', '2', '13', '5', 
                                                                          '0', '6', '1', '14', 
                                                                          '4', '11', '9'
))

Idents(T_cluster) <- 'seurat_clusters'

Idents(T_cluster)

DotPlot(T_cluster, assay = 'RNA', features = c('LTB', 'COTL1', 'CCR7', 'TCF7', 'SELL', 'IL7R', 'NEFL', 'TTC39C-AS1', 'SCGB3A1', 'GSTM3', 'FXYD2',
                                               'LINC02446', 'NELL2', 'ISG15', 'MX1', 'IFI6', 'IFIT1', 'IFI44L', 'GZMK', 'DUSP2', 'CMC1',
                                               'CCL4', 'FGFBP2', 'GZMB', 'GZMH', 'GNLY', 'NKG7', 'KLRD1', 'HOPX', 'LGALS1', 'PLEK', 'CTSW', 'TRAV17','TRBV28',
                                               'PTGDS', 'TYROBP', 'KLRC2', 'KLRC3', 'KLRB1', 'TRDC', 'TRGC1'
), 
        cols = c("white", "red")) + xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

##Dot Plot
T_cluster$sample_id <- factor(T_cluster$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                  'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                  'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                  'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(T_cluster) <- 'sample_id'


DotPlot(T_cluster, assay = 'RNA', features = c('ELK1', 'SF3B1', 'CEP55', 'NEK2', 'TOP2A', 'CHML', 'RRM2', 'PBK', 'TRIP13', 'KIF4A', 'CDC20'), 
        cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

######
# T cells_Sub
######
Sub_T_cluster_0 <- subset(T_cluster, idents = c('0'))

##Dot Plot
Sub_T_cluster_0$sample_id <- factor(Sub_T_cluster_0$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                              'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                              'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                              'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(Sub_T_cluster_0) <- 'sample_id'


DotPlot(Sub_T_cluster_0, assay = 'RNA', features = c('GZMK', 'CD69', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQB1', 'TIGIT', 'LAG3'
), 
cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

###
Sub_T_cluster_1 <- subset(T_cluster, idents = c('1'))

##Dot Plot by sample_id
Sub_T_cluster_1$sample_id <- factor(Sub_T_cluster_1$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                          'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                          'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                          'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(Sub_T_cluster_1) <- 'sample_id'


DotPlot(Sub_T_cluster_1, assay = 'RNA', features = c('GZMB', 'GZMH', 'GZMA', 'GNLY', 
                                                     'PRF1', 'NKG7', 'KLRK1', 'HLA-DPA1', 
                                                     'HLA-DPB1', 'HLA-DQB1', 'TIGIT', 'LAG3'), 
cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

###
Sub_T_cluster_2 <- subset(T_cluster, idents = c('2'))

##Dot Plot
Sub_T_cluster_2$sample_id <- factor(Sub_T_cluster_2$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                          'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                          'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                          'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(Sub_T_cluster_2) <- 'sample_id'


DotPlot(Sub_T_cluster_2, assay = 'RNA', features = c('IL2RA', 'FOXP3', 'CTLA4' ), 
cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

## markers
Idents(Final_Clus) <- 'seurat_clusters'
Idents(Final_Clus)
T_cluster <- subset(Final_Clus, idents = c('1'))
Idents(T_cluster) <- 'state'
Idents(T_cluster)

T_cluster <- JoinLayers(T_cluster, assay = 'RNA')
DEG_MGUSvsHD <- FindMarkers(T_cluster, assay = 'RNA', slot = 'data',
                            ident.1 = 'HD',
                            test.use = 'MAST',
                            min.pct = 0.2, only.pos = FALSE, verbose = TRUE)

DEG_MGUSvsHD$gene <- rownames(DEG_MGUSvsHD)

write.csv(DEG_MGUSvsHD, file = 'DEG_T_cluster_HDvsDIS.csv', row.names = FALSE) ## write CSV

Idents(T_cluster) <- 'state'
DefaultAssay(T_cluster) <- 'RNA'
VlnPlot(T_cluster, features = c("STAT1"), split.by = 'state', pt.size=0)

##Dot Plot
T_cluster$sample_id <- factor(T_cluster$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                              'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 
                                                              'MGUS_5', 'SMM_1',
                                                              'MGUS_6','MM_1',
                                                              'SMM_2', 'MM_2', 'SMM_3',
                                                              'MM_3', 'sMM_4', 'MM_4'))
Idents(T_cluster) <- 'sample_id'


DotPlot(T_cluster, assay = 'RNA', features = c('GZMH', 'PRF1', 'IFITM1', 'STAT1','IFI44L','ISG15',
                                               'MX1', 'ITGB1', 'IFI6', 'GZMK', 'GNLY',
                                               'HSPH1', 'H2AFX', 'LYZ', 'GADD45B', 
                                               'HSPA1A', 'HSPA1B'), 
cols = c("white", "red"),group.by = 'sample_id') + coord_flip() +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )



##Gene expression
DefaultAssay(sub_cluster) <- 'RNA'
FeaturePlot(sub_cluster, features = c('IFI44L'), split.by = 'state')
VlnPlot(sub_cluster, features = c("PTPRC"))
vlnplot_data <- as.data.frame(vlnplot)

##Dot Plot
DotPlot(ImG_cluster, assay = 'RNA', features = c("CD4","RPS26", "CCR7", "LTB", "TSHZ2", "ADTRP", 
                                                 "ZNF90", "MYC", "IL6ST", "MAL", "C1orf56", "HNRNPH1", 
                                                 "LRRC75A", "PHKG1", "CDC42SE1","SET", "C16orf54", "MDM4", "TNRC6B", 
                                                 "PIM2","S100A4","TNFRSF4","KLRB1","ANXA1","CRIP1","TNFAIP3", "ITGB1",
                                                 "GPR183","FOS","CRH","KRT1","RRS1","IL4R","RNF125","KPNA6","CD69","RTKN2",
                                                 "ANXA2","IL32","IFI44L","IFIT3","MX1","ISG15","XAF1","IFI44",
                                                 "IRF7","IFITM1","STAT1","CMC1","GZMK","CCL5","CST7","RGS1","CCL4","XCL2",
                                                 "DUSP2","ZFP36","CCL3","GZMH","NKG7","GNLY","FGFBP2","GZMA","GZMB","KLRD1",
                                                 "PRF1","CX3CR1","HOPX","C12orf75","TRDC","KLRC2","KLRF1","CD8A","CTSW","CD8B"), 
        cols = c("white", "red"), group.by = 'sample_id', ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )




##Gene expression
DefaultAssay(cart.combined) <- 'RNA'
FeaturePlot(cart.combined, features = c('HBB'), split.by = 'state',
            ncol = 2)
FeaturePlot(cart.combined, features = c('HBA1'))

VlnPlot(sub_cluster, features = c("NEK2"),split.by = 'state')
vlnplot_data <- as.data.frame(vlnplot)

# SaveH5Seurat(cart.combined, filename = 'results/integrated_pre_treat')

SaveH5Seurat(Final_Clus, filename = paste0(Project(object = Final_Clus), ".h5seurat"),
             overwrite = FALSE,
             verbose = TRUE,)

# top20.markers
all.markers <- FindAllMarkers(object = Final_Clus, test.use = 'MAST', min.pct = 0.1, only.pos = TRUE)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20.markers, file = 'top20.markers.csv', row.names = FALSE)

rep_cell_counts <- group_by(cart.combined@meta.data, seurat_clusters, sample_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "Cell counts.csv",sep = ",",row.names = FALSE)


##
# CD14 Mono cells
##

CD14Mono_cluster <- subset(Final_Clus, idents = c('0'))
CD14Mono_cluster <- RunPCA(CD14Mono_cluster, npcs = 50, verbose = FALSE)
CD14Mono_cluster <- FindNeighbors(CD14Mono_cluster, reduction = "pca", dims = 1:50)
CD14Mono_cluster <- FindClusters(CD14Mono_cluster, graph.name = 'integrated_snn', resolution = 0.2) # resolution: number of clusters
CD14Mono_cluster <- RunUMAP(CD14Mono_cluster, reduction = "pca", dims = 1:50, min.dist = 0.5, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(CD14Mono_cluster, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(CD14Mono_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)
DimPlot(CD14Mono_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 5)

##################
### subset to identify different group
##################

MM_cluster_subset <- subset(CD14Mono_cluster, idents = c('HD','MGUS','SMM'))

##########
# find DEG
##########

Idents(Final_Clus) <- 'seurat_clusters'
Idents(Final_Clus)
CD14Mono_cluster <- subset(Final_Clus, idents = c('0'))
Idents(CD14Mono_cluster) <- 'state'
Idents(CD14Mono_cluster)

## markers
CD14Mono_cluster <- JoinLayers(CD14Mono_cluster, assay = 'RNA')
DEG_MGUSvsHD <- FindMarkers(CD14Mono_cluster, assay = 'RNA', slot = 'data',
                            ident.1 = 'HD',
                            test.use = 'MAST',
                            min.pct = 0.2, only.pos = FALSE, verbose = TRUE)

DEG_MGUSvsHD$gene <- rownames(DEG_MGUSvsHD)

write.csv(DEG_MGUSvsHD, file = 'DEG_DISvsHD.csv', row.names = FALSE) ## write CSV

DefaultAssay(CD14Mono_cluster) <- 'RNA'
VlnPlot(CD14Mono_cluster, features = c("CD14"), split.by = 'state')


##DotPlot
CD14Mono_cluster$sample_id <- factor(CD14Mono_cluster$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                  'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                  'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                  'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(CD14Mono_cluster) <- 'sample_id'


DotPlot(CD14Mono_cluster, assay = 'RNA', features = c('IFI44L', 'IFI6', 'IFITM3', 'ISG15', 'MX1', 
                                                      'OAS3', 'IFI44', 'IRF7', 'OAS2', 'HSPB1', 
                                                      'HSP90AA1', 'HSPH1', 'DNAJB1', 'HSPA1A', 'HSPA1B'
), 
        cols = c("white", "red")) + coord_flip() +  xlab(NULL) +  ylab(NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

Idents(CD14Mono_cluster) <- 'state'
DefaultAssay(CD14Mono_cluster) <- 'RNA'
VlnPlot(CD14Mono_cluster, features = c("IFI44L"), split.by = 'state')


##
# NK cells
##

NK_cluster <- subset(Final_Clus, idents = c('4'))
NK_cluster <- RunPCA(NK_cluster, npcs = 50, verbose = FALSE)
NK_cluster <- FindNeighbors(NK_cluster, reduction = "pca", dims = 1:50)
NK_cluster <- FindClusters(NK_cluster, graph.name = 'integrated_snn', resolution = 0.2) # resolution: number of clusters
NK_cluster <- RunUMAP(NK_cluster, reduction = "pca", dims = 1:50, min.dist = 1.2, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(NK_cluster, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(NK_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)
DimPlot(NK_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 5)

DefaultAssay(NK_cluster) <- 'RNA'
FeaturePlot(NK_cluster, features = c('CD3E'))


##cell_counts

rep_cell_counts <- group_by(NK_cluster@meta.data, seurat_clusters, sample_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "NK_cluster counts.csv",sep = ",",row.names = FALSE)

# top20.markers
Idents(NK_cluster) <- 'seurat_clusters'
Idents(NK_cluster)

NK_cluster <- JoinLayers(NK_cluster, assay = 'RNA')

all.markers <- FindAllMarkers(object = NK_cluster, test.use = 'MAST', min.pct = 0.2, only.pos = TRUE)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(top20.markers, file = 'NK_cluster_top100.markers.csv', row.names = FALSE)

## markers

Idents(NK_cluster) <- 'state'
Idents(NK_cluster)

## markers
NK_cluster <- JoinLayers(NK_cluster, assay = 'RNA')
DEG_MGUSvsHD <- FindMarkers(NK_cluster, assay = 'RNA', slot = 'data',
                            ident.1 = 'HD',
                            test.use = 'MAST',
                            min.pct = 0.2, only.pos = FALSE, verbose = TRUE)

DEG_MGUSvsHD$gene <- rownames(DEG_MGUSvsHD)

write.csv(DEG_MGUSvsHD, file = 'NK_cluster_DEG_DISvsHD.csv', row.names = FALSE) ## write CSV

DefaultAssay(NK_cluster) <- 'RNA'
VlnPlot(NK_cluster, features = c("IGKC"), split.by = 'state')

##DotPlot
NK_cluster$sample_id <- factor(NK_cluster$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                          'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                          'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                          'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(NK_cluster) <- 'sample_id'


DotPlot(NK_cluster, assay = 'RNA', features = c('IFI44L', 'MX1', 'IFI6', 'ITGB7', 'ITGAX', 'IFITM1', 'IFITM3', 'ISG15', 'IRF7',
                                                'HAVCR2', 'GZMB', 'GZMH',
                                                'CX3CR1', 'CXCR4', 'CCL3', 'XCL1',
                                                'GZMK', 'LYZ',
                                                'HSPH1', 'HSPE1', 'GADD45B', 'DNAJB1', 'HSPA1A', 'HSPA1B'), 
cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

##Gene UMAP

DefaultAssay(NK_cluster) <- 'RNA'
FeaturePlot(NK_cluster, features = c('CD3G')) & theme(text = element_text(size=10),
                                                     axis.title = element_text(size=12),
                                                     legend.text=element_text(size=12),
                                                     legend.title=element_text(size=12)
)

Idents(NK_cluster) <- 'state'
DefaultAssay(NK_cluster) <- 'RNA'
VlnPlot(NK_cluster, features = c("IFI44L"), split.by = 'state')



##
# mDC cells
##

mDC_cluster <- subset(Final_Clus, idents = c('9'))
mDC_cluster <- RunPCA(mDC_cluster, npcs = 50, verbose = FALSE)
mDC_cluster <- FindNeighbors(mDC_cluster, reduction = "pca", dims = 1:50)
mDC_cluster <- FindClusters(mDC_cluster, graph.name = 'integrated_snn', resolution = 0.5) # resolution: number of clusters
mDC_cluster <- RunUMAP(mDC_cluster, reduction = "pca", dims = 1:50, min.dist = 1.2, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(mDC_cluster, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(mDC_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)
DimPlot(mDC_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 5)


# top20.markers
Idents(mDC_cluster) <- 'seurat_clusters'
Idents(mDC_cluster)

mDC_cluster <- JoinLayers(mDC_cluster, assay = 'RNA')

all.markers <- FindAllMarkers(object = mDC_cluster, test.use = 'MAST', min.pct = 0.2, only.pos = TRUE)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(top20.markers, file = 'mDC_cluster_top100.markers.csv', row.names = FALSE)


DefaultAssay(mDC_cluster) <- 'RNA'
FeaturePlot(mDC_cluster, features = c('CLEC4A')) & theme(text = element_text(size=10),
                                                     axis.title = element_text(size=12),
                                                     legend.text=element_text(size=12),
                                                     legend.title=element_text(size=12)
)

rep_cell_counts <- group_by(mDC_cluster@meta.data, seurat_clusters, sample_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "mDC_cluster counts.csv",sep = ",",row.names = FALSE)


##########
# find DEG
##########

Idents(Final_Clus) <- 'seurat_clusters'
Idents(Final_Clus)
mDC_cluster <- subset(Final_Clus, idents = c('9'))
Idents(mDC_cluster) <- 'state'
Idents(mDC_cluster)

## markers
mDC_cluster <- JoinLayers(mDC_cluster, assay = 'RNA')
DEG_MGUSvsHD <- FindMarkers(mDC_cluster, assay = 'RNA', slot = 'data',
                            ident.1 = 'HD',
                            test.use = 'MAST',
                            min.pct = 0.2, only.pos = FALSE, verbose = TRUE)

DEG_MGUSvsHD$gene <- rownames(DEG_MGUSvsHD)

write.csv(DEG_MGUSvsHD, file = 'DEG_mDC_cluster_HDvsDIS.csv', row.names = FALSE) ## write CSV

DefaultAssay(mDC_cluster) <- 'RNA'
VlnPlot(mDC_cluster, features = c("CD14"), split.by = 'state')

##DotPlot
mDC_cluster$sample_id <- factor(mDC_cluster$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                            'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                            'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                            'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(mDC_cluster) <- 'sample_id'


DotPlot(mDC_cluster, assay = 'RNA', features = c('OAS2', 'IFITM3', 'ISG15', 'IFI44', 
                                                 'IFI6', 'MX1', 'IFITM1', 'IFI44L', 
                                                 'CCL3', 'CXCL8', 'CCL3L1', 'HLA-DQA2', 
                                                 'MRC1', 'FCGR2B', 'IL1R2', 'TNFAIP3', 
                                                 'IER5', 'HSPH1', 'GADD45B', 'IER3', 
                                                 'DNAJB1', 'HSPA6', 'HSPA1A', 'HSPA1B'), 
cols = c("white", "red")) + coord_flip() + xlab(NULL) +  ylab(NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels by 45 degrees
    )



####################
#### select data
####################
Idents(sub_cluster) <- 'seurat_clusters'
Idents(sub_cluster)
sub_Tcell <- subset(sub_cluster, idents = c('3'))
DimPlot(sub_Tcell, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)

####################
#### find DEG
####################
Idents(sub_Tcell) <- 'state'
Idents(sub_Tcell)

## markers_2
markers_2 <- FindMarkers(sub_Tcell, assay = 'RNA', slot = 'data',
                         ident.1 = 'HD', ident.2 = 'MM',
                         test.use = 'MAST',
                         min.pct = 0.1, only.pos = FALSE, verbose = TRUE)
markers_2$gene <- rownames(markers_2)

write.csv(markers_2, file = 'M_2.csv', row.names = FALSE) ## write CSV

DefaultAssay(sub_Tcell) <- 'RNA'
VlnPlot(sub_Tcell, features = c("IL7R"), split.by = 'state')

#################################################
###### Subcluster analysis
#################################################

# MM cells

sub_cluster_MM <- subset(cart.combined, idents = c('10'))
sub_cluster_MM <- RunPCA(sub_cluster_MM, npcs = 50, verbose = FALSE)
sub_cluster_MM <- FindNeighbors(sub_cluster_MM, reduction = "pca", dims = 1:50)
sub_cluster_MM <- FindClusters(sub_cluster_MM, graph.name = 'integrated_snn', resolution = 0.1) # resolution: number of clusters
sub_cluster_MM <- RunUMAP(sub_cluster_MM, reduction = "pca", dims = 1:50, min.dist = 0.6, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(sub_cluster_MM, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(sub_cluster_MM, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)
DimPlot(sub_cluster_MM, reduction = "umap", label = T, label.size = 3, repel = F,
        split.by = 'sample_id', ncol = 1)


# Sub cluster_MM top50.markers
all.markers <- FindAllMarkers(object = sub_cluster_MM, test.use = 'MAST', min.pct = 0.1, only.pos = TRUE)
write.table(all.markers,file = "ABC.csv",sep = ",",row.names = FALSE)
top50.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50.markers, file = 'MM cells top50.markers.csv', row.names = FALSE)

##Sub cluster Counting
rep_cell_counts <- group_by(sub_cluster@meta.data, seurat_clusters, sample_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "sub-cluster T Cell counts.csv",sep = ",",row.names = FALSE)

##Gene expression
DefaultAssay(sub_cluster_MM) <- 'RNA'
FeaturePlot(sub_cluster_MM, features = c('NCAM1'), label = T, label.size = 3,  repel = F,
            split.by = 'sample_id', ncol = 3)
VlnPlot(sub_cluster_MM, features = c("NCAM1"))
vlnplot_data <- as.data.frame(vlnplot)

##Dot Plot
DotPlot(sub_cluster, assay = 'RNA', features = c("CD4","CD8B"), 
        cols = c("white", "red"), group.by = 'seurat_clusters', ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )


####################
#### select data
####################
Idents(sub_cluster) <- 'seurat_clusters'
Idents(sub_cluster)
sub_Tcell <- subset(sub_cluster, idents = c('3'))
DimPlot(sub_Tcell, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)

####################
#### find DEG
####################
Idents(sub_Tcell) <- 'state'
Idents(sub_Tcell)

## markers_2
markers_2 <- FindMarkers(sub_Tcell, assay = 'RNA', slot = 'data',
                         ident.1 = 'HD', ident.2 = 'MM',
                         test.use = 'MAST',
                         min.pct = 0.1, only.pos = FALSE, verbose = TRUE)
markers_2$gene <- rownames(markers_2)

write.csv(markers_2, file = 'M_2.csv', row.names = FALSE) ## write CSV

DefaultAssay(sub_Tcell) <- 'RNA'
VlnPlot(sub_Tcell, features = c("IL7R"), split.by = 'state')





#######

DefaultAssay(cart.combined) <- 'RNA'

SeuratDisk::SaveH5Seurat(cart.combined, filename = 'integrated_pre_treat_test.h5seurat')

rm(cart.list)

gc()

SaveH5Seurat(cart.combined, filename = paste0(Project(object = cart.combined), ".h5seurat"),
             overwrite = FALSE,
             verbose = TRUE,)

save(cart.combined,file="~/top.5.salaries.RData")


#######
# CD16 Mono cells
######

CD16Mono_cluster <- subset(Final_Clus, idents = c('6'))
CD16Mono_cluster <- RunPCA(CD16Mono_cluster, npcs = 50, verbose = FALSE)
CD16Mono_cluster <- FindNeighbors(CD16Mono_cluster, reduction = "pca", dims = 1:50)
CD16Mono_cluster <- FindClusters(CD16Mono_cluster, graph.name = 'integrated_snn', resolution = 0.3) # resolution: number of clusters
CD16Mono_cluster <- RunUMAP(CD16Mono_cluster, reduction = "pca", dims = 1:50, min.dist = 0.5, n.neighbors = 50) # min.dist: cluster tightness, n.neighbors: structure of cluster

DimPlot(CD16Mono_cluster, reduction = "umap",label = T, label.size = 4, repel = F)
DimPlot(CD16Mono_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)
DimPlot(CD16Mono_cluster, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'sample_id', ncol = 5)


##########
# find DEG
##########

Idents(Final_Clus) <- 'seurat_clusters'
Idents(Final_Clus)
CD16Mono_cluster <- subset(Final_Clus, idents = c('6'))
Idents(CD16Mono_cluster) <- 'state'
Idents(CD16Mono_cluster)

## markers
CD16Mono_cluster <- JoinLayers(CD16Mono_cluster, assay = 'RNA')
DEG_MGUSvsHD <- FindMarkers(CD16Mono_cluster, assay = 'RNA', slot = 'data',
                            ident.1 = 'HD',
                            test.use = 'MAST',
                            min.pct = 0.2, only.pos = FALSE, verbose = TRUE)

DEG_MGUSvsHD$gene <- rownames(DEG_MGUSvsHD)

write.csv(DEG_MGUSvsHD, file = 'CD16Mono_cluster_DEG_HDvsDIS.csv', row.names = FALSE) ## write CSV

DefaultAssay(CD14Mono_cluster) <- 'RNA'
VlnPlot(CD14Mono_cluster, features = c("CD14"), split.by = 'state')

##Dot Plot
CD16Mono_cluster$sample_id <- factor(CD16Mono_cluster$sample_id, levels = c('HD_1', 'HD_2', 'HD_3', 'HD_4', 'HD_5',
                                                                          'MGUS_1', 'MGUS_2', 'MGUS_3', 'MGUS_4', 'MGUS_5', 'MGUS_6',
                                                                          'SMM_1', 'SMM_2', 'SMM_3', 'SMM_4',
                                                                          'MM_1', 'MM_2', 'MM_3', 'MM_4'))
Idents(CD16Mono_cluster) <- 'sample_id'


DotPlot(CD16Mono_cluster, assay = 'RNA', features = c('IFITM1', 'IFIT1', 'IFI44L', 'IFI6', 'MX1', 'ISG15', 'IFI44', 'IFIT3', 'IFITM3', 'OAS3',
                                                      'OAS2', 'IRF7', 'HLA-DQA2', 'NFKBIA', 'TNF', 'IL13RA1', 'CXCL8', 'CCL3', 'IL1B', 'NLRP3',
                                                      'CCL4', 'TNFAIP3', 'CCL3L1', 'CCL4L2', 'CXCR4', 'DNAJB1', 'HSPA6', 'HSPA1A', 'HSPA1B'), 
cols = c("white", "red")) + coord_flip()  +  xlab(NULL) +  ylab(NULL) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels by 45 degrees
  )

###
# Sub cluster top20.markers
###

Idents(CD16Mono_cluster) <- 'seurat_clusters'
Idents(CD16Mono_cluster)
CD16Mono_cluster <- JoinLayers(CD16Mono_cluster, assay = 'RNA')
all.markers <- FindAllMarkers(object = CD16Mono_cluster, test.use = 'MAST', min.pct = 0.1, 
                              min.cells.feature = 100,
                              min.cells.group = 100,
                              only.pos = TRUE)

write.csv(all.markers, file = 'CD16Mono_cluster_All.csv', row.names = FALSE)


##Gene UMAP

DefaultAssay(CD16Mono_cluster) <- 'RNA'
FeaturePlot(CD16Mono_cluster, features = c('PPM1N')) & theme(text = element_text(size=10),
                                                      axis.title = element_text(size=12),
                                                      legend.text=element_text(size=12),
                                                      legend.title=element_text(size=12)
)


#######
# B cells
######

B_cluster <- subset(Final_Clus, idents = c('3'))

Idents(B_cluster) <- 'state'
DefaultAssay(B_cluster) <- 'RNA'
VlnPlot(B_cluster, features = c("IFI44L"), split.by = 'state')
