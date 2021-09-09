library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
listAttributes(ensembl)

#10x data
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

#Non-10x data
All_genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','mgi_symbol','chromosome_name','start_position','end_position'), mart = ensembl)

HugoEnsemblename=rownames(Hugosdata)
EnsembleSplit=strsplit(HugoEnsemblename,split="[.]")
Ensemblenames_withoutisoform=unlist(lapply(EnsembleSplit,'[',1))
AllHugoGenelist=match(Ensemblenames_withoutisoform,All_genes[,1])
unlabeledHugo=Ensemblenames_withoutisoform[is.na(match(Ensemblenames_withoutisoform,All_genes[,1]))]

x=All_genes[str_detect(All_genes[,3],"mir"),3]

E14Data=read.table(gzfile("P0Data.csv.gz"),header=TRUE)   
E14Data_str=CreateSeuratObject(counts=E14Data[,-1])
unlabeled=rownames(E14Data)[is.na(match(rownames(E14Data),All_genes[,3]))]

P0Data=read.table(gzfile("P0Data.csv.gz"),header=TRUE)   
P0Data_Str=CreateSeuratObject(counts=P0Data[,-1])


library(tidyr)

#E14 Data from Loo et al 2019
#P0 Data from same paper

rownames(Hugosdata)=Hugosdata[,1]

HugosDataObject=CreateSeuratObject(counts=Hugosdata[,-1])
HugosDataObject[["Percent"]]<-PercentageFeatureSet(HugosDataObject, pattern = "^MT-")

#10xdata, use pbmc instead of HugosDataObject
VlnPlot(HugosDataObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(HugosDataObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HugosDataObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
plot1a <- FeatureScatter(E14Data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2a <- FeatureScatter(E14Data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1a+plot2a
plot1b <- FeatureScatter(P0Data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2b <- FeatureScatter(P0Data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1b+plot2b



HugosDataObject_Normalized <- NormalizeData(HugosDataObject, normalization.method = "LogNormalize", scale.factor = 10000)
E14_Normalized  <- NormalizeData(E14Data, normalization.method = "LogNormalize", scale.factor = 10000)
P0_Normalized  <- NormalizeData(P0Data, normalization.method = "LogNormalize", scale.factor = 10000)

HugosDataObject_VariableFeature <- FindVariableFeatures(HugosDataObject_Normalized, selection.method = "vst", nfeatures = 2000)
E14_Feature <- FindVariableFeatures(E14_Normalized, selection.method = "vst", nfeatures = 2000)
P0_Feature <- FindVariableFeatures(P0_Normalized, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(HugosDataObject_VariableFeature), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(HugosDataObject_VariableFeature)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

E14_top10 <- head(VariableFeatures(E14_Feature), 10)
plot1 <- VariableFeaturePlot(E14_Feature)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

P0_top10 <- head(VariableFeatures(P0_Feature), 10)
plot1 <- VariableFeaturePlot(P0_Feature)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(HugosDataObject_VariableFeature)
HugosDataObject_VariableFeature_Scaled <- ScaleData(HugosDataObject_VariableFeature, features = all.genes, npcs=48 )#Default is 50 for npcs.
#The results of this are stored in pbmc[["RNA"]]@scale.data
E14.genes <- rownames(E14_Feature)
E14_Scaled <- ScaleData(E14_Feature, features = E14.genes, npcs=48 )#Default is 50 for npcs.
#The results of this are stored in pbmc[["RNA"]]@scale.data

#Dimention Reduction
HugosDataObject_VariableFeature_Scaled_PCA <- RunPCA(HugosDataObject_VariableFeature_Scaled, features = VariableFeatures(object = HugosDataObject_VariableFeature_Scaled))
E14_PCA<- RunPCA(E14_Scaled, features = VariableFeatures(object = E14_Scaled))



print(HugosDataObject_VariableFeature_Scaled_PCA[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(HugosDataObject_VariableFeature_Scaled_PCA, dims = 1:2, reduction = "pca")
DimPlot(HugosDataObject_VariableFeature_Scaled_PCA, reduction = "pca")
DimHeatmap(HugosDataObject_VariableFeature_Scaled_PCA, dims = 1, cells = 49, balanced = TRUE)
HugosDataObject_VariableFeature_Scaled_PCA_JS <- JackStraw(HugosDataObject_VariableFeature_Scaled_PCA, num.replicate = 100)
HugosDataObject_VariableFeature_Scaled_PCA_JS_Score <- ScoreJackStraw(HugosDataObject_VariableFeature_Scaled_PCA_JS, dims = 1:20)
E14_JS <- JackStraw(E14_PCA, num.replicate = 100)
E14_JS_Score <- ScoreJackStraw(E14_JS, dims = 1:20)



JackStrawPlot(HugosDataObject_VariableFeature_Scaled_PCA_JS_Score, dims = 1:15)
ElbowPlot(HugosDataObject_VariableFeature_Scaled_PCA_JS_Score)


# Cluster find neighbors

HugosDataObject_VariableFeature_Scaled_PCA_JS_Score <- FindNeighbors(HugosDataObject_VariableFeature_Scaled_PCA_JS_Score, dims = 1:10)
E14_Neighbor <- FindNeighbors(E14_JS_Score, dims = 1:10)

HugosDataObject_VariableFeature_Scaled_PCA_JS_Score <- FindClusters(HugosDataObject_VariableFeature_Scaled_PCA_JS_Score, resolution = 0.5)
E14_Cluster <- FindClusters(E14_Neighbor, resolution = 0.5)


head(Idents(HugosDataObject_VariableFeature_Scaled_PCA_JS_Score), 5)


#Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
HugosDataObject_VariableFeature_Scaled_PCA_JS_Score <- RunUMAP(HugosDataObject_VariableFeature_Scaled_PCA_JS_Score, dims = 1:10)
E14_UMAP <- RunUMAP(E14_Cluster, dims = 1:10)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(HugosDataObject_VariableFeature_Scaled_PCA_JS_Score, reduction = "umap")


#Lets use other's Data

pbmc.data <- Read10X(data.dir = "/projects/ps-zhenglab/Hugo/SequencingData/Analysis/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
