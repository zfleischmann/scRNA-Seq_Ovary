library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(devtools)
library(umap)
library(ggplot2)
library(sctransform)
library(gplots)

#Workflow is modelled after the Seurat SCTransform with integration vignette

#Large files... Otherwise integration fails
options(future.globals.maxSize=2621440000)
setwd("~/Documents/2020Wagner/CellR_Count")

#Create Objects: Read in and filter unsorted scRNA-Seq Data - code is the same until after integration
ovar.data_063 <- Read10X(data.dir="./CSEC/cellR/10X_18_PD_03_R26_063/outs/filtered_feature_bc_matrix")
ovar_CSEC <- CreateSeuratObject(counts = ovar.data_063, min.cells= 3, min.features= 200, project = "CSEC")
ovar_CSEC[["percent.mt"]] <-PercentageFeatureSet(object = ovar_CSEC, pattern = "^MT-")
ovar_CSEC <-subset(x = ovar_CSEC, subset = nFeature_RNA>200&nFeature_RNA<7000&percent.mt<25)

ovar.data_062 <- Read10X(data.dir="./GRP2/10X_18_PD_03_R26_062/outs/filtered_feature_bc_matrix/")
ovar_GRP <- CreateSeuratObject(counts = ovar.data_062, min.cells= 3, min.features= 200, project = "GRP")
ovar_GRP[["percent.mt"]] <-PercentageFeatureSet(object = ovar_GRP, pattern = "^MT-")
ovar_GRP <-subset(x = ovar_GRP, subset = nFeature_RNA>200&nFeature_RNA<7000&percent.mt<25)


### SCtransform integration and normalization:
#Pre-processing for integration - SCTransform will issue 50+ warnings. These can be (safely) ignored per advice from the developers.
ovar_CSEC <- PercentageFeatureSet(ovar_CSEC, pattern = "^MT-", col.name = "percent.mt")
ovar_CSEC <- SCTransform(ovar_CSEC, vars.to.regress = "percent.mt", verbose = FALSE)
ovar_GRP <- PercentageFeatureSet(ovar_GRP, pattern = "^MT-", col.name = "percent.mt")
ovar_GRP <- SCTransform(ovar_GRP, vars.to.regress = "percent.mt", verbose = FALSE)

#Integrate samples
Ovar.list <- c(ovar_CSEC, ovar_GRP)
Ovar.features <- SelectIntegrationFeatures(object.list=Ovar.list, nfeatures=3000)
Ovar.list <- PrepSCTIntegration(object.list = Ovar.list, anchor.features=Ovar.features, verbose=FALSE)
Ovar_unsorted.anchors <- FindIntegrationAnchors(object.list = Ovar.list, normalization.method = "SCT", anchor.features = Ovar.features, verbose=FALSE)
Ovar_unsorted.integrated <- IntegrateData(anchorset = Ovar_unsorted.anchors, normalization.method="SCT", verbose=FALSE)
DefaultAssay(object = Ovar_unsorted.integrated)

### Original Normalization 
#Processing for integration: Original normalization from Wagner 2020
#ovar_unsorted_CSEC <- NormalizeData(object = ovar_CSEC, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#ovar_unsorted_CSEC <- FindVariableFeatures(object = ovar_unsorted_CSEC, selection.method = "vst",nfeatures = 2000, verbose = FALSE)
#ovar_unsorted_GRP <- NormalizeData(object = ovar_GRP, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#ovar_unsorted_GRP <- FindVariableFeatures(object = ovar_unsorted_GRP, selection.method = "vst",nfeatures = 2000, verbose = FALSE)
#DefaultAssay(object=Ovar_unsorted.integrated) <- 'integrated'

########################################################################################################################
# Dimensionality reduction and clustering
# Per the vignette: DO NOT run ScaleData after integration of sctransform normalized seurat objects. 
# Instead, proceed with dimensionality reduction and clustering.
########################################################################################################################
Ovar_unsorted.integrated <- RunPCA(Ovar_unsorted.integrated, verbose = FALSE)
Ovar_unsorted.integrated <- RunUMAP(Ovar_unsorted.integrated, dims = 1:30, verbose = FALSE)
Ovar_unsorted.integrated <- FindNeighbors(Ovar_unsorted.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
Ovar_unsorted.integrated <- FindClusters(Ovar_unsorted.integrated, verbose = FALSE, resolution = 0.1)

#Plot of clusters
Ovar_unsorted.integrated <- RenameIdents(object = Ovar_unsorted.integrated, "0"="1", "1"="2", "2"="3", "3"="4", "4"="5", "5"="6", "6"="7", "7"="8", "8"="9", "9"="10")
d <- DimPlot(Ovar_unsorted.integrated, label = TRUE, repel=TRUE) + NoLegend()
  
########################################################################################################################
# Gene Expression Scores (see Wagner et al 2020)
########################################################################################################################

DefaultAssay(object = Ovar_unsorted.integrated) <- "SCT"
stroma <- c("DCN", 'PDGFRA', 'APOE', 'FHL2')
PV <- c('RGS5', 'MCAM', 'RERGL', 'TAGLN', 'MYH11')
ENDO <- c('CD34', 'VWF', 'FLI1', 'CDH5')
Gran <- c('AMH', 'FST', 'FOXL2', 'BEX1')
IMMUNE <- c('CD69', 'ITGB2', 'CXCR4', 'CD14')
oo <- c("GDF9", 'ZP3', 'OOSP2', 'FIGLA')
osc <- c("DAZL", "DPPA3", "PRDM1")

oo_percent <- 100*colSums(Ovar_unsorted.integrated@assays$RNA@counts[oo,])/colSums(Ovar_unsorted.integrated@assays$RNA@counts)
Ovar_unsorted.integrated <- AddMetaData(object = Ovar_unsorted.integrated, metadata = oo_percent, col.name='Oocyte_percent')
#oo <- VlnPlot(Ovar_unsorted.integrated, features=c("Oocyte_percent"), pt.size=0.05) + NoLegend() + stat_summary(fun = median, geom='point', size = 2, colour='black')

osc <- c("DAZL", "DPPA3", "PRDM1")
osc_percent <- 100*colSums(Ovar_unsorted.integrated@assays$RNA@counts[osc,])/colSums(Ovar_unsorted.integrated@assays$RNA@counts)
Ovar_unsorted.integrated <- AddMetaData(object = Ovar_unsorted.integrated, metadata = osc_percent, col.name='OSC_percent')
#osc <- VlnPlot(Ovar_unsorted.integrated, features=c("OSC_percent"), pt.size=0.05) + NoLegend() + stat_summary(fun = median, geom='point', size = 2, colour='black')

oocyte_subset <- subset(Ovar_unsorted.integrated, idents = 10)

f <- FeatureScatter(oocyte_subset, feature1="OSC_percent", feature2="Oocyte_percent", pt.size = 2)
f <- f + labs( title = "",
               x = "OSC expression: DAZL, DPPA3, PRDM1",
               y = "Oocyte expression: GDF9, ZP3, OOSP2, FIGLA")+ NoLegend()

########################################################################################################################
# Combine plots: Figure 1A+B
########################################################################################################################

d+f