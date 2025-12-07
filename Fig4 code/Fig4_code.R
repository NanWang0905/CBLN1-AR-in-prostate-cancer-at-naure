## Set system language to English
Sys.setenv(LANGUAGE = "en")

## Disable automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear working environment
rm(list=ls())

setwd("F:/mydata/GSE210358/")
getwd()

library(Seurat)
library(tidyverse)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(harmony)
library(DoubletFinder)


#### Data loading ####
# Sample names and file paths
samples <- list(
  HMP04 = "rawdata_human/GSM6428952_1778_JZ_HMP_04_IGO_10726_2_dense.csv.gz",
  HMP05 = "rawdata_human/GSM6428953_1779_HMP05_IGO_10726_3_dense.csv.gz",
  HMP08 = "rawdata_human/GSM6428954_1845_HMP-08_IGO_10837_19_dense.csv.gz",
  HMP11_1 = "rawdata_human/GSM6428955_1968_HMP11_1_IGO_11247_3_dense.csv.gz",
  HMP11_2 = "rawdata_human/GSM6428956_1969_HMP11_2_IGO_11247_4_dense.csv.gz",
  HMP13 = "rawdata_human/GSM6428957_1846_JZHP_3_IGO_10837_21_dense.csv.gz",
  HMP14 = "rawdata_human/GSM6428958_2284_HMP_14B_IGO_11588_19_dense.csv.gz",
  HMP15 = "rawdata_human/GSM6428959_2331_HMP_15_IGO_11377_B_22_dense.csv.gz",
  HMP16 = "rawdata_human/GSM6428960_2016_JZHP_04_IGO_11245_7_dense.csv.gz",
  HMP17 = "rawdata_human/GSM6428961_2513_HMP_17_IGO_11874_35_dense.csv.gz",
  HMP19 = "rawdata_human/GSM6428962_2518_HMP_19_IGO_11874_45_dense.csv.gz",
  HMP20 = "rawdata_human/GSM6428963_2624_HMP20_IGO_12065_37_dense.csv.gz",
  HMP25 = "rawdata_human/GSM6428964_3226_SZ-1327_HMP25_IGO_12437_254_dense.csv.gz",
  HMP26 = "rawdata_human/GSM6428965_3309_SZ-1357_HMP26_IGO_12437_I_19_dense.csv.gz"
)

# Loop through each sample and create Seurat objects
for (name in names(samples)) {
  file <- samples[[name]]
  
  # Read
  tmp <- t(read.csv(file, row.names = 1))
  colnames(tmp) <- paste0(name, "_", 1:ncol(tmp))
  tmp <- as.data.frame(tmp)
  tmp <- tmp[-1, ]
  
  # Create Seurat object
  seu <- CreateSeuratObject(tmp, project = name)
  
  # Save object under its sample name
  assign(name, seu)
}

#### Doublet removal ####
# Remove doublets

# HMP04
HMP04 <- NormalizeData(HMP04)
HMP04 <- FindVariableFeatures(HMP04, selection.method = "vst", nfeatures = 2000)
HMP04 <- ScaleData(HMP04)
HMP04 <- RunPCA(HMP04)
HMP04 <- RunUMAP(HMP04, dims = 1:10)
# Estimated doublet rate ~2%
nExp_poi <- round(0.02*nrow(HMP04@meta.data))
# Identify doublets
HMP04 <- doubletFinder_v3(HMP04, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP04$DF.classifications_0.25_0.09_82)
DimPlot(HMP04,group.by = "DF.classifications_0.25_0.09_82")
# Remove doublets using subset
HMP04=subset(HMP04, subset = DF.classifications_0.25_0.09_82 != "Doublet")

# HMP05
HMP05 <- NormalizeData(HMP05)
HMP05 <- FindVariableFeatures(HMP05, selection.method = "vst", nfeatures = 2000)
HMP05 <- ScaleData(HMP05)
HMP05 <- RunPCA(HMP05)
HMP05 <- RunUMAP(HMP05, dims = 1:10)
# Estimated doublet rate ~2.7%
nExp_poi <- round(0.027*nrow(HMP05@meta.data))
# Identify doublets
HMP05 <- doubletFinder_v3(HMP05, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP05$DF.classifications_0.25_0.09_151)
DimPlot(HMP05,group.by = "DF.classifications_0.25_0.09_151")
# Remove doublets
HMP05=subset(HMP05, subset = DF.classifications_0.25_0.09_151 != "Doublet")

# HMP08
HMP08 <- NormalizeData(HMP08)
HMP08 <- FindVariableFeatures(HMP08, selection.method = "vst", nfeatures = 2000)
HMP08 <- ScaleData(HMP08)
HMP08 <- RunPCA(HMP08)
HMP08 <- RunUMAP(HMP08, dims = 1:10)
# Estimated doublet rate ~3.9%
nExp_poi <- round(0.039*nrow(HMP08@meta.data))
# Identify doublets
HMP08 <- doubletFinder_v3(HMP08, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP08$DF.classifications_0.25_0.09_315)
DimPlot(HMP08,group.by = "DF.classifications_0.25_0.09_315")
# Remove doublets
HMP08=subset(HMP08, subset = DF.classifications_0.25_0.09_315 != "Doublet")

# HMP11_1
HMP11_1 <- NormalizeData(HMP11_1)
HMP11_1 <- FindVariableFeatures(HMP11_1, selection.method = "vst", nfeatures = 2000)
HMP11_1 <- ScaleData(HMP11_1)
HMP11_1 <- RunPCA(HMP11_1)
HMP11_1 <- RunUMAP(HMP11_1, dims = 1:10)
# Estimated doublet rate ~3.8%
nExp_poi <- round(0.038*nrow(HMP11_1@meta.data))
# Identify doublets
HMP11_1 <- doubletFinder_v3(HMP11_1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP11_1$DF.classifications_0.25_0.09_293)
DimPlot(HMP11_1,group.by = "DF.classifications_0.25_0.09_293")
# Remove doublets
HMP11_1=subset(HMP11_1, subset = DF.classifications_0.25_0.09_293 != "Doublet")

# HMP11_2
HMP11_2 <- NormalizeData(HMP11_2)
HMP11_2 <- FindVariableFeatures(HMP11_2, selection.method = "vst", nfeatures = 2000)
HMP11_2 <- ScaleData(HMP11_2)
HMP11_2 <- RunPCA(HMP11_2)
HMP11_2 <- RunUMAP(HMP11_2, dims = 1:10)
# Estimated doublet rate ~4.3%
nExp_poi <- round(0.043*nrow(HMP11_2@meta.data))
# Identify doublets
HMP11_2 <- doubletFinder_v3(HMP11_2, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP11_2$DF.classifications_0.25_0.09_385)
DimPlot(HMP11_2,group.by = "DF.classifications_0.25_0.09_385")
# Remove doublets
HMP11_2=subset(HMP11_2, subset = DF.classifications_0.25_0.09_385 != "Doublet")

# HMP13
HMP13 <- NormalizeData(HMP13)
HMP13 <- FindVariableFeatures(HMP13, selection.method = "vst", nfeatures = 2000)
HMP13 <- ScaleData(HMP13)
HMP13 <- RunPCA(HMP13)
HMP13 <- RunUMAP(HMP13, dims = 1:10)
# Estimated doublet rate ~1.9%
nExp_poi <- round(0.019*nrow(HMP13@meta.data))
# Identify doublets
HMP13 <- doubletFinder_v3(HMP13, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP13$DF.classifications_0.25_0.09_74)
DimPlot(HMP13,group.by = "DF.classifications_0.25_0.09_74")
# Remove doublets
HMP13=subset(HMP13, subset = DF.classifications_0.25_0.09_74 != "Doublet")

# HMP14
HMP14 <- NormalizeData(HMP14)
HMP14 <- FindVariableFeatures(HMP14, selection.method = "vst", nfeatures = 2000)
HMP14 <- ScaleData(HMP14)
HMP14 <- RunPCA(HMP14)
HMP14 <- RunUMAP(HMP14, dims = 1:10)
# Estimated doublet rate ~2.7%
nExp_poi <- round(0.027*nrow(HMP14@meta.data))
# Identify doublets
HMP14 <- doubletFinder_v3(HMP14, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP14$DF.classifications_0.25_0.09_152)
DimPlot(HMP14,group.by = "DF.classifications_0.25_0.09_152")
# Remove doublets
HMP14=subset(HMP14, subset = DF.classifications_0.25_0.09_152 != "Doublet")

# HMP15
HMP15 <- NormalizeData(HMP15)
HMP15 <- FindVariableFeatures(HMP15, selection.method = "vst", nfeatures = 2000)
HMP15 <- ScaleData(HMP15)
HMP15 <- RunPCA(HMP15)
HMP15 <- RunUMAP(HMP15, dims = 1:10)
# Estimated doublet rate ~1.2%
nExp_poi <- round(0.012*nrow(HMP15@meta.data))
# Identify doublets
HMP15 <- doubletFinder_v3(HMP15, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP15$DF.classifications_0.25_0.09_32)
DimPlot(HMP15,group.by = "DF.classifications_0.25_0.09_32")
# Remove doublets
HMP15=subset(HMP15, subset = DF.classifications_0.25_0.09_32 != "Doublet")

# HMP16
HMP16 <- NormalizeData(HMP16)
HMP16 <- FindVariableFeatures(HMP16, selection.method = "vst", nfeatures = 2000)
HMP16 <- ScaleData(HMP16)
HMP16 <- RunPCA(HMP16)
HMP16 <- RunUMAP(HMP16, dims = 1:10)
# Estimated doublet rate ~3.1%
nExp_poi <- round(0.031*nrow(HMP16@meta.data))
# Identify doublets
HMP16 <- doubletFinder_v3(HMP16, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP16$DF.classifications_0.25_0.09_197)
DimPlot(HMP16,group.by = "DF.classifications_0.25_0.09_197")
# Remove doublets
HMP16=subset(HMP16, subset = DF.classifications_0.25_0.09_197 != "Doublet")

# HMP17
HMP17 <- NormalizeData(HMP17)
HMP17 <- FindVariableFeatures(HMP17, selection.method = "vst", nfeatures = 2000)
HMP17 <- ScaleData(HMP17)
HMP17 <- RunPCA(HMP17)
HMP17 <- RunUMAP(HMP17, dims = 1:10)
# Estimated doublet rate ~2.6%
nExp_poi <- round(0.026*nrow(HMP17@meta.data))
# Identify doublets
HMP17 <- doubletFinder_v3(HMP17, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP17$DF.classifications_0.25_0.09_139)
DimPlot(HMP17,group.by = "DF.classifications_0.25_0.09_139")
# Remove doublets
HMP17=subset(HMP17, subset = DF.classifications_0.25_0.09_139 != "Doublet")

# HMP19
HMP19 <- NormalizeData(HMP19)
HMP19 <- FindVariableFeatures(HMP19, selection.method = "vst", nfeatures = 2000)
HMP19 <- ScaleData(HMP19)
HMP19 <- RunPCA(HMP19)
HMP19 <- RunUMAP(HMP19, dims = 1:10)
# Estimated doublet rate ~3.1%
nExp_poi <- round(0.031*nrow(HMP19@meta.data))
# Identify doublets
HMP19 <- doubletFinder_v3(HMP19, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP19$DF.classifications_0.25_0.09_195)
DimPlot(HMP19,group.by = "DF.classifications_0.25_0.09_195")
# Remove doublets
HMP19=subset(HMP19, subset = DF.classifications_0.25_0.09_195 != "Doublet")

# HMP20
HMP20 <- NormalizeData(HMP20)
HMP20 <- FindVariableFeatures(HMP20, selection.method = "vst", nfeatures = 2000)
HMP20 <- ScaleData(HMP20)
HMP20 <- RunPCA(HMP20)
HMP20 <- RunUMAP(HMP20, dims = 1:10)
# Estimated doublet rate ~0.3%
nExp_poi <- round(0.003*nrow(HMP20@meta.data))
# Identify doublets
HMP20 <- doubletFinder_v3(HMP20, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Examine classification
table(HMP20$DF.classifications_0.25_0.09_2)
DimPlot(HMP20,group.by = "DF.classifications_0.25_0.09_2")
# Remove doublets
HMP20=subset(HMP20, subset = DF.classifications_0.25_0.09_2 != "Doublet")

# HMP25
HMP25 <- NormalizeData(HMP25)
HMP25 <- FindVariableFeatures(HMP25, selection.method = "vst", nfeatures = 2000)
HMP25 <- ScaleData(HMP25)
HMP25 <- RunPCA(HMP25)
HMP25 <- RunUMAP(HMP25, dims = 1:10)
# Approximately 2%
nExp_poi <- round(0.02*nrow(HMP25@meta.data))
# Identify doublets
HMP25 <- doubletFinder_v3(HMP25, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Check doublet classification
table(HMP25$DF.classifications_0.25_0.09_83)
DimPlot(HMP25,group.by = "DF.classifications_0.25_0.09_83")
# Remove doublets using subset
HMP25 <- subset(HMP25, subset = DF.classifications_0.25_0.09_83 != "Doublet")


# HMP26
HMP26 <- NormalizeData(HMP26)
HMP26 <- FindVariableFeatures(HMP26, selection.method = "vst", nfeatures = 2000)
HMP26 <- ScaleData(HMP26)
HMP26 <- RunPCA(HMP26)
HMP26 <- RunUMAP(HMP26, dims = 1:10)
# Approximately 1.4%
nExp_poi <- round(0.014*nrow(HMP26@meta.data))
# Identify doublets
HMP26 <- doubletFinder_v3(HMP26, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# Check doublet classification
table(HMP26$DF.classifications_0.25_0.09_40)
DimPlot(HMP26,group.by = "DF.classifications_0.25_0.09_40")
# Remove doublets
HMP26 <- subset(HMP26, subset = DF.classifications_0.25_0.09_40 != "Doublet")


#### merge & QC ####

# Merge samples
mCRPC <- merge(HMP04, y = c(HMP05,HMP08,HMP11_1,HMP11_2,HMP13,HMP14,HMP15,HMP16,
                            HMP17,HMP19,HMP20,HMP25,HMP26),
               project = "mCRPC")
table(mCRPC$orig.ident)

# Remove raw objects
rm(HMP04,HMP05,HMP08,HMP11_1,HMP11_2,HMP13,HMP14,HMP15,HMP16,HMP17,HMP19,HMP20,HMP25,HMP26)

# Calculate mitochondrial percentage
mCRPC[["percent.mt"]] <- PercentageFeatureSet(mCRPC, pattern = "^MT.")

# Calculate ribosomal percentage
mCRPC[["percent.rb"]] <- PercentageFeatureSet(mCRPC, pattern = c("^RP[SL]"))

# Check basic metrics
VlnPlot(mCRPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 3)

# QC filtering
mCRPC.qc <- subset(mCRPC, subset = nFeature_RNA > 500 & percent.mt < 20 & nCount_RNA > 1000)
rm(mCRPC)
table(mCRPC.qc@meta.data$orig.ident)


#### Basic processing ####
# Normalization, HVG selection, scaling, PCA
mCRPC.qc <- NormalizeData(mCRPC.qc, verbose = TRUE)
mCRPC.qc <- FindVariableFeatures(mCRPC.qc, selection.method = "vst", nfeatures = 2000)
mCRPC.qc <- ScaleData(mCRPC.qc, verbose = TRUE, vars.to.regress = "percent.rb")  # Regress out ribosomal genes
mCRPC.qc <- RunPCA(mCRPC.qc, npcs = 50, verbose = TRUE)
# Batch correction using Harmony
mCRPC.qc <- RunHarmony(mCRPC.qc, group.by.vars = 'orig.ident', reduction = "pca", reduction.save = "harmony")
# Clustering
mCRPC.qc <- FindNeighbors(mCRPC.qc, dims = 1:15, reduction = "harmony")
mCRPC.qc <- FindClusters(mCRPC.qc, resolution = 0.5)  # Higher resolution → more clusters
# Dimensionality reduction
mCRPC.qc <- RunUMAP(mCRPC.qc, dims = 1:20, reduction = "harmony")
mCRPC.qc <- RunTSNE(mCRPC.qc, dims = 1:20, reduction = "harmony")
# Visualization
DimPlot(mCRPC.qc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(mCRPC.qc, reduction = "umap", group.by = "orig.ident")


#### Annotation ####
# Check markers
DotPlot(mCRPC.qc, features = c(
  "EPCAM","KRT8","KRT5","CDH1",      # epithelial
  "ENG","CLDN5","VWF","CDH5",        # endothelial
  "DCN","TNFAIP6","APOD","FBLN1",    # fibroblasts
  "CD79A","CD79B","CD19","MS4A1",    # B cells
  "CD2","CD3D","CD3E","CD3G","NCAM1","FCGR3A",  # T/NK
  "CD14","CD68","AIF1","CSF1R",      # myeloid
  "MS4A2","ENPP3","FCER1A","KIT"     # mast cells
))

# Assign cell types
mCRPC.qc@meta.data[which(mCRPC.qc@meta.data$seurat_clusters %in% c(0,11,23,25)), "celltype"] <- "T/NK cells"
mCRPC.qc@meta.data[which(mCRPC.qc@meta.data$seurat_clusters %in% c(2,21,24)), "celltype"] <- "Myeloid"
mCRPC.qc@meta.data[which(mCRPC.qc@meta.data$seurat_clusters %in% c(5,8,22)), "celltype"] <- "Fibroblasts"
mCRPC.qc@meta.data[which(mCRPC.qc@meta.data$seurat_clusters %in% c(6,26)), "celltype"] <- "Endothelial"
mCRPC.qc@meta.data[which(mCRPC.qc@meta.data$seurat_clusters %in% c(18)), "celltype"] <- "B cells"
mCRPC.qc@meta.data[which(mCRPC.qc@meta.data$seurat_clusters %in% c(20)), "celltype"] <- "Mast"
mCRPC.qc@meta.data[which(mCRPC.qc@meta.data$seurat_clusters %in% c(1,3,4,7,9,10,12:17,19)), "celltype"] <- "Epithelial"

# Label metastatic sites
table(mCRPC.qc$orig.ident)
mCRPC.qc@meta.data[mCRPC.qc$orig.ident %in% c("HMP04","HMP05","HMP16","HMP25","HMP26"),
                   "tissue"] <- "liver"
mCRPC.qc@meta.data[mCRPC.qc$orig.ident %in% c("HMP13","HMP14","HMP19","HMP20"),
                   "tissue"] <- "LN"
mCRPC.qc@meta.data[mCRPC.qc$orig.ident %in% c("HMP11_1","HMP11_2"),
                   "tissue"] <- "brain"
mCRPC.qc@meta.data[mCRPC.qc$orig.ident %in% c("HMP15"),
                   "tissue"] <- "Epidural"
mCRPC.qc@meta.data[mCRPC.qc$orig.ident %in% c("HMP08"),
                   "tissue"] <- "Mediastinal mass"
mCRPC.qc@meta.data[mCRPC.qc$orig.ident %in% c("HMP17"),
                   "tissue"] <- "Left retrocrural abdomen mass"
table(mCRPC.qc$tissue)

# Label NEPC
mCRPC.qc@meta.data[, "PCa_type"] <- "CRPC"
mCRPC.qc@meta.data[mCRPC.qc$orig.ident %in% c("HMP04","HMP16","HMP17"),
                   "PCa_type"] <- "NEPC"

table(mCRPC.qc$orig.ident, mCRPC.qc$PCa_type)


#### Remove NEPC subclusters ####

# Subset
non_NEPC = mCRPC.qc[, mCRPC.qc@meta.data$PCa_type %in% c("CRPC")]

# Processing
non_NEPC <- NormalizeData(non_NEPC, verbose = TRUE)
non_NEPC <- FindVariableFeatures(non_NEPC, selection.method = "vst", nfeatures = 2000)
non_NEPC <- ScaleData(non_NEPC, verbose = TRUE, vars.to.regress = "percent.rb")
non_NEPC <- RunPCA(non_NEPC, npcs = 50, verbose = TRUE)
non_NEPC <- RunHarmony(non_NEPC, group.by.vars = 'orig.ident', reduction = "pca", reduction.save = "harmony")
non_NEPC <- FindNeighbors(non_NEPC, dims = 1:15, reduction = "harmony")
non_NEPC <- FindClusters(non_NEPC, resolution = 0.5)
non_NEPC <- RunUMAP(non_NEPC, dims = 1:15, reduction = "harmony")

# Visualization
DimPlot(non_NEPC, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(non_NEPC, reduction = "umap", group.by = 'orig.ident')

# Marker checking
DotPlot(mCRPC.qc, features = c(
  "EPCAM","KRT8","KRT5","CDH1",
  "ENG","CLDN5","VWF","CDH5",
  "DCN","TNFAIP6","APOD","FBLN1",
  "CD79A","CD79B","CD19","MS4A1",
  "CD2","CD3D","CD3E","CD3G","NCAM1","FCGR3A",
  "CD14","CD68","AIF1","CSF1R",
  "MS4A2","ENPP3","FCER1A","KIT"
))

# Re-annotation
non_NEPC@meta.data[which(non_NEPC@meta.data$seurat_clusters %in% c(1,9,13)), "celltype"] <- "T/NK"
non_NEPC@meta.data[which(non_NEPC@meta.data$seurat_clusters %in% c(2,21)), "celltype"] <- "myeloid"
non_NEPC@meta.data[which(non_NEPC@meta.data$seurat_clusters %in% c(5,8)), "celltype"] <- "fibroblast"
non_NEPC@meta.data[which(non_NEPC@meta.data$seurat_clusters %in% c(7)), "celltype"] <- "endo"
non_NEPC@meta.data[which(non_NEPC@meta.data$seurat_clusters %in% c(14)), "celltype"] <- "Bcell"
non_NEPC@meta.data[which(non_NEPC@meta.data$seurat_clusters %in% c(17)), "celltype"] <- "mast"
non_NEPC@meta.data[which(non_NEPC@meta.data$seurat_clusters %in% c(0,3,4,6,10,11,12,15,16,18,19,20,22)), "celltype"] <- "epithelial"

# Annotation
DimPlot(non_NEPC, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 7) +
  ggtitle("CellType") +
  theme(plot.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 15))

# Gene expression visualization
FeaturePlot(non_NEPC, features = "NRXN1", cols = c("grey", "red"), pt.size = 1, order = TRUE)
FeaturePlot(non_NEPC, features = "AR", cols = c("grey", "red"), pt.size = 1)
FeaturePlot(non_NEPC, features = "KLK3", cols = c("grey", "red"), pt.size = 1)
FeaturePlot(non_NEPC, features = "FKBP5", cols = c("grey", "red"), pt.size = 1)
FeaturePlot(non_NEPC, features = "TMPRSS2", cols = c("grey", "red"), pt.size = 1)

# Save data
saveRDS(non_NEPC, file = "non_NEPC.rds")
