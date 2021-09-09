BiocManager::install("devtools")
BiocManager::install("dendextend")
BiocManager::install("dplyr")
BiocManager::install("matrixStats")
BiocManager::install("Matrix")
#BiocManager::install("limma")
BiocManager::install("AnnotationDbi")
BiocManager::install("GO.db")
BiocManager::install("preprocessCore")
BiocManager::install("impute")
BiocManager::install("DESeq2")
devtools::install_github("AllenInstitute/tasic2016data")
devtools::install_github("AllenInstitute/scrattch.hicat")

library(limma)
library(tasic2016data)
library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

str(tasic_2016_anno)

head(tasic_2016_anno)
dim(tasic_2016_counts)

tasic_2016_counts[1:5,1:5]

select.cells <- tasic_2016_anno %>%
  filter(primary_type_label != "unclassified") %>%
  filter(grepl("Igtp|Ndnf|Vip|Sncg|Smad3",primary_type_label)) %>%
  select(sample_name) %>%
  unlist()

length(select.cells)

ref_anno <- tasic_2016_anno %>%
  filter(sample_name %in% select.cells)

# Make a data.frame of unique cluster id, type, color, and broad type
ref.cl.df <- ref_anno %>%
  select(primary_type_id, 
         primary_type_label, 
         primary_type_color, 
         broad_type) %>%
  unique()

#standardize cluster annoation with cluster_id, cluster_label and cluster_color. These are the required fields to visualize clusters properly.
colnames(ref.cl.df)[1:3] <- c("cluster_id", "cluster_label", "cluster_color")

# Sort by cluster_id
ref.cl.df <- arrange(ref.cl.df, cluster_id)
row.names(ref.cl.df) <- ref.cl.df$cluster_id

# Set up the ref.cl factor object
ref.cl <- setNames(factor(ref_anno$primary_type_id), ref_anno$sample_id)

#Covert to CPM
tasic_2016_cpm <- cpm(tasic_2016_counts[,select.cells])
norm.dat <- log2(tasic_2016_cpm + 1)
norm.dat <- Matrix(norm.dat, sparse = TRUE)

#Set parameter
de.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.5, #0.3 if low depth
                     q2.th       = NULL,
                     q.diff.th   = 1, #splitting cell types based on graded or combinatorial differences 0 
                     de.score.th = 40, #Large dataset(>10000), 150
                     min.cells = 4)
#Dimention filter
gene.counts <- colSums(norm.dat > 0)
rm.eigen <- matrix(log2(gene.counts), ncol = 1)
row.names(rm.eigen) <- names(gene.counts)
colnames(rm.eigen) <- "log2GeneCounts"

#clustering
strict.param <- de_param(de.score.th = 500)

onestep.result <- onestep_clust(norm.dat, 
                                select.cells = select.cells, 
                                dim.method = "pca", #or WGCNA
                                de.param = strict.param, 
                                rm.eigen = rm.eigen)

pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
display.result <- display_cl(onestep.result$cl, norm.dat, plot = TRUE, de.param = de.param)
dev.off()

#Iteration
iter.result <- iter_clust(norm.dat, 
                          select.cells = select.cells, 
                          dim.method = "pca", 
                          de.param = de.param, 
                          rm.eigen = rm.eigen,
                          result = onestep.result)

#Test
dim.method <- "WGCNA"

explore.param <- de_param(padj.th     = 0.05, 
                          lfc.th      = 1, 
                          low.th      = 1, 
                          q1.th       = 0.5,
                          q2.th       = NULL,
                          q.diff.th   = 0.7, 
                          de.score.th = 40)

explore.result <- onestep_clust(norm.dat, 
                                select.cells = select.cells, 
                                dim.method = dim.method, 
                                de.param = explore.param, 
                                rm.eigen = rm.eigen)

display.result <- display_cl(explore.result$cl, 
                             norm.dat, 
                             plot = TRUE, 
                             de.param = explore.param)


## Merging

rd.dat <- t(norm.dat[iter.result$markers, select.cells])

merge.param <- de_param(de.score.th = 70) # The original value was 40.

merge.result <- merge_cl(norm.dat, 
                         cl = iter.result$cl, 
                         rd.dat = rd.dat,
                         de.param = merge.param)

display.result <- display_cl(merge.result$cl, 
                             norm.dat, 
                             plot = TRUE, 
                             de.param = merge.param)

# Set up the cl and cl.df objects for use with compare_annotate()
iter.cl <- setNames(as.factor(iter.result$cl), select.cells)
iter.cl.df <- data.frame(cluster_id = unique(iter.cl),
                         cluster_label = paste0("Pre-merge_cl_",unique(iter.cl)),
                         cluster_color = rainbow(length(unique(iter.cl))))
rownames(iter.cl.df) <- iter.cl.df$cluster_id

compare.result <- compare_annotate(merge.result$cl, iter.cl, iter.cl.df)
compare.result$g

# Generate comparison
compare.result <- compare_annotate(iter.result$cl, ref.cl, ref.cl.df)
# Output the plot
compare.result$g

# Get cl factors and data.frame 
cl <- compare.result$cl
cl.df <- compare.result$cl.df
display.result <- display_cl(cl, 
                             norm.dat, 
                             plot=TRUE, 
                             de.param = de.param,
                             n.markers = 20)

de.genes <- display.result$de.genes

set.seed(12345)
result <- run_consensus_clust(norm.dat[,select.cells], 
                              select.cells = select.cells,
                              niter = 100, 
                              de.param = de.param, 
                              rm.eigen = rm.eigen, 
                              dim.method = "pca", 
                              output_dir = "subsample_PCA")
#Mydata

setwd("/projects/ps-zhenglab/Hugo/SequencingData/HugosData/Analysis/Regen_NoRegen1/")
Hugosdata=read.table("Combined_RNAseq.txt",header=TRUE)
Hugosdata_cpm=cpm(as.matrix(Hugosdata[,-1]))
names(Hugosdata_cpm)=Hugosdata[,1]
Hnorm.dat <- log2(Hugosdata_cpm + 1)
Hnorm.dat <- Matrix(Hnorm.dat, sparse = TRUE)
#Set parameter
Hde.param <- de_param(padj.th     = 0.05, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.5, #0.3 if low depth
                     q2.th       = NULL,
                     q.diff.th   = 1, #splitting cell types based on graded or combinatorial differences 0 
                     de.score.th = 40, #Large dataset(>10000), 150
                     min.cells = 4)
#Dimention filter

Hgene.counts <- colSums(Hnorm.dat > 0)
Hrm.eigen <- matrix(log2(Hgene.counts), ncol = 1)
row.names(Hrm.eigen) <- names(gene.counts)
colnames(Hrm.eigen) <- "log2GeneCounts"
#clustering
Hstrict.param <- de_param(de.score.th = 500)

Honestep.result <- onestep_clust(Hnorm.dat, 
                                select.cells = Hselect.cells, 
                                dim.method = "pca", #or WGCNA
                                de.param = Hstrict.param, 
                                rm.eigen = Hrm.eigen)

pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
display.result <- display_cl(Honestep.result$cl, norm.dat, plot = TRUE, de.param = de.param)
dev.off()

#Iteration
iter.result <- iter_clust(Hnorm.dat, 
                          select.cells = Hselect.cells, 
                          dim.method = "pca", 
                          de.param = Hde.param, 
                          rm.eigen = Hrm.eigen,
                          result = onestep.result)

