########  Single-nuclei Multiome (GEX + ATAC data) Analysis of Canine OSA {Kyland} ########  
### Written By: Rebecca L. Nance, 04/20/2023 ###

#adapted from the following vignettes: 
  #Seurat Guided Clustering Tutorial (GEX only): https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
  #Seurat Integrating scRNAseq and scATACseq data: https://satijalab.org/seurat/articles/atacseq_integration_vignette.html 
  #Seurat Weighted Nearest Neighbor Analysis: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html 
  #Signac Joint RNA and ATAC analysis: https://stuartlab.org/signac/articles/pbmc_multiomic.html

###### 1. Prepare Workspace (seurat) ############ 
# load libraries and functions
library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)
library(HGNChelper)
library(openxlsx)
library(ggplot2)
library(Signac)
library(data.table)
library(magrittr)
library(ggraph)
library(igraph)
library(tidyverse)
library(data.tree)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggplot2)
setwd("C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/CellRangerOutput")
#load the data (this folder contains the 'filtered_feature_bc_matrix.h5' file)
counts <- Read10X_h5('C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/CellRangerOutput/filtered_feature_bc_matrix.h5')
fragpath <- "C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/CellRangerOutput/atac_fragments.tsv.gz"
#create a Seurat object containing the GEX(RNA) data
sn_OSA <- CreateSeuratObject(counts = counts$'Gene Expression', assay="RNA", project="sn_OSA")
#ensure the GEX(RNA) data is loaded
Assays(sn_OSA)
#get summary stats on the RNA(GEX) data
ngenes=colSums(GetAssayData(object=sn_OSA, assay="RNA",slot="counts")>0)
summary(ngenes)
mean(ngenes)

#create ATAC assay and add it to the Seurat object
sn_OSA[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath)
#ensure the RNAseq data and ATACdata is loaded in the Seurat assays
Assays(sn_OSA)
#get summary stats on the ATAC data
ngenes=colSums(GetAssayData(object=sn_OSA, assay="ATAC", slot="counts")>0)
summary(ngenes)
mean(ngenes)


###### 2. Quality Control and Pre-Processing (seurat) ############ 
#FIGURE2A: VIOLIN PLOT QC METRICS--Visualize QC metrics as a violin plot
p1 <- VlnPlot(sn_OSA, features = c("nFeature_RNA"), pt.size = 0)
p1 <- p1 + labs(title="GEX: Features/Nucleus") + 
  scale_y_continuous(labels=scientific, n.breaks=4) +
  theme(axis.title.x=element_blank(), plot.title=element_text(size=12) )
p2 <- VlnPlot(sn_OSA, features = c("nCount_RNA"), pt.size = 0)
p2 <- p2 + labs(title="GEX: Counts/Nucleus") + 
  scale_y_continuous(labels=scientific, n.breaks=4) +
  theme(axis.title.x=element_blank(), plot.title=element_text(size=12) )
p3 <- VlnPlot(sn_OSA, features = c("nFeature_ATAC"), pt.size = 0)
p3 <- p3 + labs(title="ATAC: Features/Nucleus") + 
  scale_y_continuous(labels=scientific, n.breaks=4) +
  theme(axis.title.x=element_blank(), plot.title=element_text(size=12) )
p4 <- VlnPlot(sn_OSA, features = c("nCount_ATAC"), pt.size = 0)
p4 <- p4 + labs(title="ATAC: Counts/Nucleus") + 
  scale_y_continuous(labels=scientific, n.breaks=4) +
  theme(axis.title.x=element_blank(), plot.title=element_text(size=12) )
p5 <- ggarrange(p1, p2, ncol=2, legend=FALSE)
p6 <- ggarrange(p3, p4, ncol=2, legend=FALSE)
ggarrange(p5 + p6)

#filter out nuclei based on counts (UMIs) and features (genes)
sn_OSA_filter <- subset(x=sn_OSA, 
                        subset = nFeature_RNA > 100 & nFeature_RNA < 30000 &
                          nCount_RNA > 50 & nCount_RNA < 50000 &
                          nFeature_ATAC > 100 & nFeature_ATAC < 30000 &
                          nCount_ATAC > 50 & nCount_ATAC < 50000)
#get summary stats on the filtered RNA(GEX) data
ngenes=colSums(GetAssayData(object=sn_OSA_filter, assay="RNA", slot="counts")>0)
summary(ngenes)
mean(ngenes)
#get summary stats on the filtered ATAC data
ngenes=colSums(GetAssayData(object=sn_OSA_filter, assay="ATAC", slot="counts")>0)
summary(ngenes)
mean(ngenes)


### Fig2B: SCATTERPLOTS--Before & After Filtering
#PLOT: before filtering--RNA counts vs RNA features
p1<-FeatureScatter(sn_OSA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ scale_x_continuous(labels=scientific, n.breaks=4) + theme(axis.text=element_text(size=10))
p1 <- p1 + labs(title="GEX Before Filtering: 0.87", x="nCount_GEX", y="nFeature_GEX")+
  theme(plot.title=element_text(size=12))
#PLOT: filtered RNA counts vs RNA features 
p2<-FeatureScatter(sn_OSA_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ scale_x_continuous(labels=scientific, n.breaks=4) + theme(axis.text=element_text(size=10))
p2 <- p2 + labs(title="GEX After Filtering: 0.95", x="nCount_GEX", y="nFeature_GEX")+
  theme(plot.title=element_text(size=12))

#PLOT: before filtering--ATAC counts vs ATAC features
p3<-FeatureScatter(sn_OSA, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")+ scale_x_continuous(labels=scientific, n.breaks=4) + theme(axis.text=element_text(size=10))
p3 <- p3+ labs(title="ATAC Before Filtering: 0.98")+
  theme(plot.title=element_text(size=12))

#PLOT: filtered ATAC counts vs ATAC features
p4<-FeatureScatter(sn_OSA_filter, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")+ scale_x_continuous(labels=scientific, n.breaks=4) + theme(axis.text=element_text(size=10))
p4 <- p4+ labs(title="ATAC After Filtering: 1")+
  theme(plot.title=element_text(size=12))
p5 <- ggarrange(p1, p2, ncol=2, legend=FALSE)
p5
p6 <- ggarrange(p3, p4, ncol=2, legend=FALSE)
p6
ggarrange(p5, p6)


###### 3. Clustering/WNN (seurat) ############ 
#GEX data
DefaultAssay(sn_OSA_filter) <- "RNA"
#sn_OSA_filter <- SCTransform(sn_OSA_filter, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
sn_OSA_filter <- SCTransform(sn_OSA_filter)
#identify the 20 most highly variable genes
top20 <- head(VariableFeatures(sn_OSA_filter), 20)
top20
#FIGURE5a: volcano plot - variable features with labels
p1 <- VariableFeaturePlot(sn_OSA_filter)
p2 <- LabelPoints(plot = p1, points = top20, repel = TRUE)
plot5a <- p2 + theme(legend.position="right")
plot5a
#run PCA and cluster the GEX data
sn_OSA_filter <-  RunPCA(sn_OSA_filter) %>% RunUMAP(dims = 1:20, reduction = 'pca', reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#ATAC data
DefaultAssay(sn_OSA_filter) <- "ATAC"
sn_OSA_filter <- RunTFIDF(sn_OSA_filter)
sn_OSA_filter <- FindTopFeatures(sn_OSA_filter, min.cutoff = 'q0')
sn_OSA_filter <- RunSVD(sn_OSA_filter)
sn_OSA_filter <- RunUMAP(sn_OSA_filter, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#Combined WNN
sn_OSA_filter <- FindMultiModalNeighbors(sn_OSA_filter, reduction.list = list("pca", "lsi"), dims.list = list(1:20, 2:20))
sn_OSA_filter <- RunUMAP(sn_OSA_filter, nn.name = "weighted.nn", reduction.name = "UMAP", reduction.key = "wnnUMAP_")
#identify 9 clusters (c0-8)
sn_OSA_filter <- FindClusters(sn_OSA_filter, graph.name = "wsnn", algorithm = 3, resolution=0.1, verbose = FALSE)

#plot all 3 UMAP graphs
plot3a <- DimPlot(sn_OSA_filter, reduction = "umap.rna", label = TRUE, repel = TRUE, label.size=6) + ggtitle("GEX")+ theme(plot.title=element_text(hjust=0.5)) + NoLegend()
plot3b <- DimPlot(sn_OSA_filter, reduction = "umap.atac", label = TRUE, repel = TRUE, label.size=6) + ggtitle("ATAC")+ theme(plot.title=element_text(hjust=0.5)) + NoLegend()
plot3c <- DimPlot(sn_OSA_filter, reduction = "UMAP", label = TRUE, repel = TRUE, label.size=6) + ggtitle("WNN")+ theme(plot.title=element_text(hjust=0.5)) 
plot3a
plot3b
plot3c
#extract the metadata
md <- sn_OSA_filter@meta.data %>% as.data.table
###GET THE NUMBER OF CELLS FROM EACH CLUSTER
table(sn_OSA_filter@meta.data$wsnn_res.0.1, sn_OSA_filter@meta.data$orig.ident)
#extract the barcodes and cluster ID [[to annotate copykat/infercnv results with]]
export_df <- sn_OSA_filter@meta.data %>% 
  rownames_to_column("barcodes") %>%
  select(barcodes, seurat_clusters)
write.table(export_df, "C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/barcodes_clusters.txt", 
            sep="\t", 
            row.names=FALSE, 
            col.names=FALSE, 
            quote = FALSE)


###### 4. Cluster Annotation with Bone/OSA/Immune Markers (scType) ###### 
#ScType for Single Cell Annotation: https://github.com/IanevskiAleksandr/sc-type 
# load gene set preparation function and cell type annotation function 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load the annotation file (DB file)
db_ = 'C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/ScType_bone_NEW.xlsx'
# identify the tissue type
tissue = "Bone" 
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = sn_OSA_filter[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(sn_OSA_filter@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sn_OSA_filter@meta.data[sn_OSA_filter@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sn_OSA_filter@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

#overlay identified cell types on UMAP plot
sn_OSA_filter@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sn_OSA_filter@meta.data$customclassif[sn_OSA_filter@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])}
ccolss=brewer.pal(n=8, name="Dark2")
#FIGURE6A--annotation with known cell markers
plot6a <- DimPlot(sn_OSA_filter, reduction = "UMAP", label = TRUE, repel = TRUE, cols = ccolss, group.by = 'customclassif') + 
  ggtitle("Cluster Annotation with Known Cell Markers") + theme(plot.title = element_text(hjust = 0, size=11))
plot6a 


###### 5. Cluster Annotation with Bulk Tumor/Normal Markers (scType) ###### 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = 'C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/sc_markers_bulkRNAseq_ScType.xlsx'
tissue = "Osteosarcoma" 
gs_list = gene_sets_prepare(db_, tissue)
es.max = sctype_score(scRNAseqData = sn_OSA_filter[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive) 
cL_resutls = do.call("rbind", lapply(unique(sn_OSA_filter@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sn_OSA_filter@meta.data[sn_OSA_filter@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sn_OSA_filter@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
sn_OSA_filter@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sn_OSA_filter@meta.data$customclassif[sn_OSA_filter@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])}
ccolss=c('#0096FF', '#FF0000', '#808080')
#FIGURE5B--bulk OSA tumor + normal bone RNAseq annotation
plot5b <- DimPlot(sn_OSA_filter, reduction = "UMAP", label = FALSE, cols = ccolss, group.by = 'customclassif') + 
  ggtitle("Cluster Annotation with Bulk OSA Tumor/Normal Bone RNAseq Markers") + 
  theme(plot.title = element_text(hjust = 0.5, size=11))
plot5b <- LabelClusters(plot5b, id="customclassif", fontface="bold", color="black", repel=TRUE)
plot5b



###### 6. Plot Specific Marker Genes (seurat) ######
DefaultAssay(sn_OSA_filter) <- "SCT"

#plot marker genes for each cluster 
#osteoblasts (clusters 0,1,7)
FeaturePlot(sn_OSA_filter, 
            features=c("RUNX2","CDH11","PCNA","ACAN","MKI67","TOP2A","COL1A1"), 
            reduction='UMAP', max.cutoff=3, ncol=4)
#fibroblasts (cluster 2)
FeaturePlot(sn_OSA_filter, 
            features=c("LUM","DCN","VIM","THY1","FAP","PRRX1","COL1A1"), 
            reduction='UMAP', max.cutoff=3, ncol=4)
#endothelial (cluster 3)
FeaturePlot(sn_OSA_filter, 
            features=c("CDH5","PECAM1","EGFL7","CD93","ENG","EMCN"),
            reduction='UMAP', max.cutoff=3, ncol=4)
#myeloid cells (cluster 4)
FeaturePlot(sn_OSA_filter, 
            features=c("CD14","CD74"), 
            reduction='UMAP', max.cutoff=3, ncol=4)
#osteoclasts (cluster 5)
FeaturePlot(sn_OSA_filter, 
            features=c("ATP6V0D2","DCSTAMP","CTSK","OCSTAMP","MMP9","ACP5"), 
            reduction='UMAP', max.cutoff=3, ncol=4)
#osteocytes (cluster 6)
FeaturePlot(sn_OSA_filter, 
            features=c("BGLAP","SPP1","CD86","IBSP"), 
            reduction='UMAP', max.cutoff=3, ncol=4)
#memory CD4+ T cells (cluster 8)
FeaturePlot(sn_OSA_filter, 
            features=c("CD4","CD3E","CD3D","CTLA4","LCK","LTB","CD2"),
            reduction='UMAP', max.cutoff=3, ncol=4)

#plot marker genes based on tumor bulk annotation data
FeaturePlot(sn_OSA_filter, 
            features=c("HOXC10","SPAG5","TOP2A","IQGAP3","HELLS","MKI67","CLSPN","RAD54L"),
            reduction='UMAP', max.cutoff=3, ncol=4)

#comparison to top genes from bulk tumor/normal RNA seq 
#top 8 upreg genes in bulk
FeaturePlot(sn_OSA_filter, features=c("GTSE1","HELLS","SPAG5","RAD54L","IQGAP3","CIT","HOXC10","TOP2A","MKI67"), reduction='UMAP', max.cutoff=3)
#top 8 upreg genes in patient C
FeaturePlot(sn_OSA_filter, features=c("TFPI2","DDX60","OAS1","CD5L","TERT","OAS2","RPGRIP1L","OAS3","ANLN"), reduction='UMAP', max.cutoff=3)



#plot M1 macrophage marker genes 
FeaturePlot(sn_OSA_filter,
            features=c("CD80","CD86","CD68"),
            reduction='UMAP', max.cutoff=3, ncol=3)
#plot M2 macrophage marker genes 
FeaturePlot(sn_OSA_filter,
            features=c("MRC1","CD163"),
            reduction='UMAP', max.cutoff=3, ncol=3)

#plot markers of CD8+ cytotoxic T cells 
FeaturePlot(sn_OSA_filter,
            features=c("CD3E","CD8A","CD8B","CD5","CD27","CD28"),
            reduction='UMAP', max.cutoff=3, ncol=3)


#plot markers of CD8+ T-cell exhaustion 
FeaturePlot(sn_OSA_filter,
            features=c("HAVCR2","LAG3","BATF","VHL","FOXO1","FOXP1"),
            reduction='UMAP', max.cutoff=3, ncol=3)


#plot marker genes based on immunotherapies
FeaturePlot(sn_OSA_filter,
            features=c("CD8A","LOX", "CD8B","NKG7","TNFRSF6"),
            reduction='UMAP', max.cutoff=3, ncol=4)

###plotting 'random' genes 
#chemo resistance gene: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2889491/
FeaturePlot(sn_OSA_filter, features=c("GSTP1"), reduction='UMAP', max.cutoff=3)
#survivin gene (associated with mets): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2889491/
FeaturePlot(sn_OSA_filter, features=c("BIRC5"), reduction='UMAP', max.cutoff=3)
#FAS gene: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2889491/
FeaturePlot(sn_OSA_filter, features=c("FAS"), reduction='UMAP', max.cutoff=3)
#osteoblast differentiation markers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3571113/#:~:text=The%20most%20frequently%20used%20markers,and%20PPR%20(Table%202).
FeaturePlot(sn_OSA_filter, features=c("BGLAP","ALPL","SPP1","IBSP"), reduction='UMAP', max.cutoff=3)


###### 7. Identification of Differentially Expressed Genes (seurat) ############ 
DefaultAssay(sn_OSA_filter) <- "SCT"
#find markers for every cluster compared to all remaining cells and report only the positive ones 
sn_OSA_filter.markers <- FindAllMarkers(sn_OSA_filter, assay="SCT", test.use="roc",  min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
#save this output to a file 
write.csv(sn_OSA_filter.markers, "C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/DifferentialAnalysisOutput/cluster_markers_Seurat_ROC_onlypos_9clusters.csv", row.names=TRUE)

#FIG5B:heatmap of the top 5 DEGs defining each cluster (based on GEX data, which was 'SCTrasform'-ed)
sn_OSA_filter.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
plot5b <- DoHeatmap(sn_OSA_filter, features = top5$gene) + theme(text=element_text(size=9), legend.position="bottom")
plot5b


### Elucidation of the Osteoblasts (DEGs between clusters 0,1,7)
#find cluster 0 DEGs distinguishing cluster 0 (OB, unknown) from clusters 1 and 7 (OB, tumor) 
cluster0vs1and7.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1 = 0, ident.2 = c(1, 7), min.pct = 0.25, logfc.threshold = 0.25)
#find clusters 1&7 DEGs distinguishing clusters 1 (OB, tumor) and 7 (tumor) from cluster 0 (OB, unknown) 
cluster1and7vs0.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1 = c(1, 7), ident.2 = 0, min.pct = 0.25, logfc.threshold = 0.25)

#find cluster 0 DEGs distinguishing cluster 0 from cluster 1
cluster0vs1.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1=0, ident.2=1, min.pct = 0.25, logfc.threshold = 0.25)
#find cluster 0 DEGs distinguishing cluster 0 from cluster 7
cluster0vs7.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1=0, ident.2=7, min.pct = 0.25, logfc.threshold = 0.25)
#find cluster 1 DEGs distinguishing cluster 1 from cluster 0
cluster1vs0.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1=1, ident.2=0, min.pct = 0.25, logfc.threshold = 0.25)
#find cluster 1 DEGs distinguishing cluster 1 from cluster 7
cluster1vs7.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1=1, ident.2=7, min.pct = 0.25, logfc.threshold = 0.25)
#find cluster 7 DEGs distinguishing cluster 7 from cluster 0
cluster7vs0.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1=7, ident.2=0, min.pct = 0.25, logfc.threshold = 0.25)
#find cluster 7 DEGs distinguishing cluster 7 from cluster 1
cluster7vs1.markers <- FindMarkers(sn_OSA_filter, assay="SCT", test.use="roc", ident.1=1, ident.2=0, min.pct = 0.25, logfc.threshold = 0.25)

#write the files
write.csv(cluster1and7v0.markers, "C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/DifferentialAnalysisOutput/cluster_0_vs_1_and_7_markers.csv")


###### 8. Panther (GO Pathway) Analysis Plots ######
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

#import the PANTHER pathway results {for each cluster}
setwd("C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/DifferentialAnalysisOutput/PANTHER_Pathway_Analysis_9clusters")
data <- read.table("cluster0_GOpathways.txt", sep="\t", fill=TRUE, header=TRUE)
#data <- read.table("cluster0vs1and7_GOpathways.txt", sep="\t", fill=TRUE, header=TRUE)
#data <- read.table("cluster1and7vs0_GOpathways.txt", sep="\t", fill=TRUE, header=TRUE)

#extract the GO terms + raw p-value + FDR + Fold enrichment
data <- data[,c(1,6,7,8)]
#extract the top 5 pathways
data <- head(data, 5)
#extract the fold enrichment 
fold <- data$Client.Text.Box.Input..fold.Enrichment.
fold <- as.numeric(fold)
#lock in the order of the GO terms based on enrichment value
data$GO.biological.process.complete <- factor(data$GO.biological.process.complete, levels=data$GO.biological.process.complete[order(data$Client.Text.Box.Input..fold.Enrichment.)])
#extract the -log10 of the p-value
counts <- -log10(data$Client.Text.Box.Input..FDR.)


#for clusters 0-9: bubble plot with legend
c0 <- ggplot(data, 
             aes(x=GO.biological.process.complete,  y=counts, fill=GO.biological.process.complete, color=GO.biological.process.complete, size=fold)) + 
  geom_point(alpha=1) + 
  scale_size(trans="log10", name="Fold Enrichment", range=c(3,4.5)) +
  coord_flip() + 
  ylab("-log10(FDR)") +
  xlab("GO Biological Process")

c0 <- c0 + ggtitle("Enriched GO Biological Processes in Cluster 0") + 
  theme(plot.title=element_text(hjust=1, size=10), axis.title=element_text(size=7), axis.text.x=element_text(size=7),axis.text.y=element_text(size=8), legend.title=element_text(size=8)) + 
  guides(color="none", fill=FALSE) 
c0

ggsave("cluster0_GOpathways.png",
       path="C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/Figures",
       width=9, height=1.5, dpi=300)


#for clusters 0 vs 1and7 and vice versa: bubble plot with legend 
c0 <- ggplot(data, 
             aes(x=GO.biological.process.complete,  y=counts, fill=GO.biological.process.complete, color=GO.biological.process.complete, size=fold)) + 
  geom_point(alpha=1) + 
  scale_size(trans="log10", name="Fold Enrichment", range=c(3,4.5)) +
  coord_flip() + 
  ylab("-log10(FDR)") +
  xlab("GO Biological Process")
c0 + ggtitle("Enriched GO Biological Processes in Cluster 0 Compared to Clusters 1 and 7") + 
  theme(plot.title=element_text(hjust=1, size=10), axis.title=element_text(size=7), axis.text.x=element_text(size=7),axis.text.y=element_text(size=8), legend.title=element_text(size=8)) + 
  guides(color="none", fill=FALSE)

#export using 650x160
#or export to create a png file: 
ggsave("cluster0_GOpathways.png",
       path="C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/Figures",
       width=8, height=2, dpi=300)


###### 9. Single Cell GSEA--Hallmark & Canonical Pathway Analysis (singleseqgset) ######
#https://arc85.github.io/singleseqgset/articles/singleseqgset.html
library(singleseqgset)
library(msigdbr)
library(heatmap3)

### Hallmark Pathways 
### download gene sets of interest using msigdbr[[DOG, Hallmark pathways]]
c.dog <- msigdbr(species="Canis lupus familiaris", category="H")
c.names <- unique(c.dog$gs_name)
c.sets <- vector("list",length=length(c.names))
names(c.sets) <- c.names
for (i in names(c.sets)) {
  c.sets[[i]] <- pull(c.dog[c.dog$gs_name==i,"gene_symbol"])
}
### extract the seurat data 
DefaultAssay(sn_OSA_filter) <- "SCT"
exp.mat <- GetAssayData(object=sn_OSA_filter)
### use singleseqgset to perform gene set enrichment analysis 
logfc.data <- logFC(cluster.ids=sn_OSA_filter@meta.data$seurat_clusters, expr.mat=exp.mat)
gse.res <- wmw_gsea(expr.mat=exp.mat,
                    cluster.cells=logfc.data[[1]],
                    log.fc.cluster=logfc.data[[2]],
                    gene.sets=c.sets)
res.stats <- gse.res[["GSEA_statistics"]]
res.pvals <- gse.res[["GSEA_p_values"]]
#correct for multiple comparisons
res.pvals <- apply(res.pvals,2,p.adjust,method="fdr")
#write the results to a file
write.csv(res.stats, "C:/users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/gsea.hallmark.stats.csv") 
write.csv(res.pvals, "C:/users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/gsea.hallmark.pvals.csv")
#top 25 gene sets enriched by Z scores (for cluster 0)
#res.stats <- res.stats[order(res.stats[,1],decreasing=TRUE)[1:25],]
#top 25 gene sets by pvalues 
#res.pvals <- res.pvals[order(res.stats[,1],decreasing=TRUE)[1:25],] 
#Fig 7A--create heatmap of the Z scores for Hallmark Pathways 
#heatmap3(res.stats, Colv=NA, cexRow=0.8, scale="row", main="Hallmark Pathways",
#         useRaster=TRUE, showRowDendro = FALSE, showColDendro = TRUE)
heatmap3(res.stats, cexRow=0.8, scale="row", main="Hallmark Pathways",
         useRaster=TRUE, showRowDendro = FALSE)

### Canonical pathways (C2=CP)
### download gene sets of interest using msigdbr[[DOG, cell type signature gene sets]]
c.dog <- msigdbr(species="Canis lupus familiaris", category="C2", subcategory="CP")
c.names <- unique(c.dog$gs_name)
c.sets <- vector("list",length=length(c.names))
names(c.sets) <- c.names
for (i in names(c.sets)) {
  c.sets[[i]] <- pull(c.dog[c.dog$gs_name==i,"gene_symbol"])
}
DefaultAssay(sn_OSA_filter) <- "SCT"
exp.mat <- GetAssayData(object=sn_OSA_filter)
### use singleseqgset to perform gene set enrichment analysis 
logfc.data <- logFC(cluster.ids=sn_OSA_filter@meta.data$seurat_clusters, expr.mat=exp.mat)
gse.res <- wmw_gsea(expr.mat=exp.mat,
                    cluster.cells=logfc.data[[1]],
                    log.fc.cluster=logfc.data[[2]],
                    gene.sets=c.sets)
res.stats <- gse.res[["GSEA_statistics"]]
res.pvals <- gse.res[["GSEA_p_values"]]
#correct for multiple comparisons
res.pvals <- apply(res.pvals,2,p.adjust,method="fdr")
#write the results to a file 
write.csv(res.stats, "C:/users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/gsea.canonical.stats.csv") 
write.csv(res.pvals, "C:/users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/gsea.canonical.pvals.csv")
#top 20 gene sets enriched by Z scores
#res.stats <- res.stats[order(res.stats[,1],decreasing=TRUE)[1:20],]
#top 20 gene sets by pvalues 
#res.pvals <- res.pvals[order(res.stats[,1],decreasing=TRUE)[1:20],] 
#Fig 7A--create heatmap of the Z scores for Hallmark Pathways 
heatmap3(res.stats, cexRow=0.8, scale="row", main="Canonical Pathways",
         useRaster=TRUE, showRowDendro = FALSE)


###### 10. Copy Number Variation Analysis  (infercnv) ######
#https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
#https://bioconductor.org/packages/release/bioc/vignettes/infercnv/inst/doc/inferCNV.html
library(infercnv)
library(stringr)
library(tidyverse)
setwd("C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB")

## 1: create gene_order_file with 4 columns ['gene symbol, chr, location start, location stop'], tab-delim, no col header
gene_ordering_file <- read.table("C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/Ref.filtered.gtf", sep="\t", fill=TRUE, header=FALSE)
#first extract only the genes
gene_ordering_file <- gene_ordering_file[str_detect(gene_ordering_file$V3, "gene"),]
#extract the last column
gene_ordering_file <- separate(data=gene_ordering_file, col=V9, into=c("a", "b"), sep=";")
gene_ordering_file <- separate(data=gene_ordering_file, col=a, into=c("a", "b", "c", "d"))
#extract only the gene name, gene ID, and start & stop locations in the chromosome 
#gene_ordering_file <- gene_ordering_file[c(11,1,4:5)]
gene_ordering_file <- gene_ordering_file[c(1,11,4:5)]
#colnames(gene_ordering_file) <- c("gene","NC_chr","start","stop")
colnames(gene_ordering_file) <- c("NC_chr","gene","start","stop")
#get annotations for the chromosomes and genbank refseq 
annotations <- read.table("C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/chr_genbankrefseq_anno.txt")
annotations <- annotations[,c(1,4)]
colnames(annotations) <- c("chr", "NC_chr")

#merge the annotations with the results
gene_ordering_file <- merge(gene_ordering_file, annotations, by="NC_chr")
gene_ordering_file <- gene_ordering_file[,c(2,5,3,4)]
#gene_ordering_file <- gene_ordering_file[,c(1,5,3,4)]
#remove duplicate rows
###THIS IS THE ISSUE  : 
#gene_ordering_file1 <- gene_ordering_file %>% distinct(c, .keep_all=TRUE) 
non_dupl_idx <- which(duplicated(gene_ordering_file$gene) == FALSE)
#non_dupl_idx <- which(duplicated(gene_ordering_file$c) == FALSE)
gene_ordering_file <- gene_ordering_file[non_dupl_idx, ]


#write it to a file
write.table(gene_ordering_file, 
            "C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/gene_ordering_file_infercnv.txt", 
            sep="\t", 
            col.names=FALSE, 
            row.names=FALSE, 
            quote=FALSE)
## 2: create sample annotation file with 2 columns['cell name(barcode), cell type(cluster ID)'], tab-delim, no col header
#extract the barcodes and new cluster ID
DefaultAssay(sn_OSA_filter) <- "SCT"
export_df <- sn_OSA_filter@meta.data %>% 
  rownames_to_column("barcodes") %>%
  select(barcodes, seurat_clusters )
write.table(export_df, "C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/annotation_file_infercnv.txt", 
            sep="\t", 
            row.names=FALSE, 
            col.names=FALSE, 
            quote = FALSE)
## 3: create the matrix file variable 
counts_matrix = as.matrix(sn_OSA_filter@assays$RNA@counts[,colnames(sn_OSA_filter)]) 

## create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file="annotation_file_infercnv.txt",
                                    gene_order_file="gene_ordering_file_infercnv.txt", 
                                    ref_group_names=c("2","3","4","5","6","8"))
## run the full default analysis 
#out_dir = "C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/infercnv_output_final"
out_dir = "C:/Users/becca/Desktop/multiome/infercnv_output_final"
infercnv_obj_default = infercnv::run(infercnv_obj,
                                     out_dir=out_dir, 
                                     cutoff=0.1, 
                                     cluster_by_groups=FALSE,
                                     denoise=TRUE,
                                     HMM=TRUE,
                                     png_res=300
)



infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff=0.1, 
  out_dir=out_dir,
  cluster_by_groups=FALSE, 
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=TRUE,
  no_prelim_plot=TRUE,
  png_res=300,
)

plot_cnv(infercnv_obj_default,
         outdir=outdir,
         title="infercnv figure",
         obs_title="Malignant cells",
         ref_title="Normal cells",
         cluster_by_groups=FALSE,
         plot_chr_scale=TRUE,
         color_safe_pal=TRUE,
         output_filename="infercnv_scaled_to_chr")


###### 11. (not using for pub) Immune Landscape: Subset the Immune Cell Clusters (clusters 4, 5, and 8:myeloid,osteoclast,& memCD4+ T cells) ###### 
subset <- subset(x=sn_OSA_filter, idents=c("4","5","8"))
##perform PCA and UMAP projections for this subset of clusters with GEX data
DefaultAssay(subset) <- "RNA"
#GEX data
DefaultAssay(subset) <- "RNA"
subset <- SCTransform(subset)
subset <-  RunPCA(subset) %>% RunUMAP(dims = 1:10, reduction = 'pca', reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#ATAC data
DefaultAssay(subset) <- "ATAC"
subset <- RunTFIDF(subset)
subset <- FindTopFeatures(subset, min.cutoff = 'q0')
subset <- RunSVD(subset)
subset <- RunUMAP(subset, reduction = 'lsi', dims = 2:10, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
#Combined WNN
subset <- FindMultiModalNeighbors(subset, reduction.list = list("pca", "lsi"), dims.list = list(1:10, 2:10))
subset <- RunUMAP(subset, nn.name = "weighted.nn", reduction.name = "UMAP", reduction.key = "wnnUMAP_")
#identify 9 clusters (c0-8)
subset <- FindClusters(subset, graph.name = "wsnn", algorithm = 3, resolution=0.7, verbose = FALSE)
#plot all 3 UMAP graphs
plot3a <- DimPlot(subset, reduction = "umap.rna", label = TRUE, repel = TRUE, label.size=6) + ggtitle("GEX")+ theme(plot.title=element_text(hjust=0.5)) + NoLegend()
plot3b <- DimPlot(subset, reduction = "umap.atac", label = TRUE, repel = TRUE, label.size=6) + ggtitle("ATAC")+ theme(plot.title=element_text(hjust=0.5)) + NoLegend()
plot3c <- DimPlot(subset, reduction = "UMAP", label = TRUE, repel = TRUE, label.size=6) + ggtitle("WNN")+ theme(plot.title=element_text(hjust=0.5)) 
plot3a
plot3b
plot3c

#ScType for Single Cell Annotation: https://github.com/IanevskiAleksandr/sc-type 
# load gene set preparation function and cell type annotation function 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load the annotation file (DB file)
db_ = 'C:/Users/becca/Dropbox/SingleNucleiMultiomeSeq/OSA_snMultiomeSeq_GITHUB/ScType_bone_NEW.xlsx'
# identify the tissue type
tissue = "Bone" 
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = subset[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(subset@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(subset@meta.data[subset@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(subset@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
#overlay identified cell types on UMAP plot
subset@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  subset@meta.data$customclassif[subset@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])}
ccolss=brewer.pal(n=8, name="Dark2")
#FIGURE6A--annotation with known cell markers
plot12a <- DimPlot(subset, reduction = "UMAP", label = TRUE, repel = TRUE, cols = ccolss, group.by = 'customclassif') + 
  ggtitle("Annotation of the Subclustered Immune Cells with Known Markers") + theme(plot.title = element_text(hjust = 0, size=11))
plot12a 

#plot specific genes(not using for pub)
DefaultAssay(subset) <- "SCT"
#plot M1 macrophage marker genes 
FeaturePlot(subset,
            features=c("CD80","CD86","CD68"),
            reduction='UMAP', max.cutoff=3, ncol=3)
#plot M2 macrophage marker genes 
FeaturePlot(subset,
            features=c("MRC1","CD163"),
            reduction='UMAP', max.cutoff=3, ncol=3)
#plot CD8 T cell marker genes
FeaturePlot(subset,
            features=c("CD8A","CD8B"),
            reduction='UMAP', max.cutoff=3, ncol=3)
