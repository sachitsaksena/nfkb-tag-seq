.#####
# EXECUTE
#####
library(dplyr)
library(tibble)
library(RColorBrewer)
library(gplots)
library(genefilter)
library(calibrate)
library(DESeq2)
library(data.table)
library(wesanderson)
# source original function list 
# source("./helper_functions.R")
source("silhouette.R")
source("rld_pca.R")
# setwd("~/2018_projects/yankeelov/tag-based_RNAseq/")

##########
# DATA PRE-PROCESSING
##########

# load original data frame
a_data <- read.csv("./A.csv", header=TRUE, row.names = 1)
countdata <- read.csv("./allcounts.csv", header=TRUE, row.names=1)

colnames(countdata) <- gsub("\\.fastq.trim.sam.counts", "", colnames(countdata))
colnames(countdata) <- gsub("X", "", colnames(countdata))
colnames(a_data) <- gsub("\\.fastq.trim.sam.counts", "", colnames(a_data))
colnames(a_data) <- gsub("X", "", colnames(a_data))

countdata_table <- data.table(countdata, keep.rownames = TRUE)
a_table <- data.table(a_data, keep.rownames = TRUE)
allcounts <- merge(countdata_table, a_table, all.x=TRUE)

# separate conditions
allcounts %>% 
  dplyr::select(rn, starts_with("2"), starts_with("A")) -> condition2
# separate conditions
allcounts %>%
  dplyr::select(rn, matches("^[A-U]")) -> condition1

colnames(condition1)[1] <- "a_0gene"
colnames(condition2)[1] <- "a_0gene"
# sort column names by letter
condition1_frame <- data.frame(condition1)
condition1_sorted <- condition1_frame[ , order(names(condition1_frame))]
# sort column names by letter 
condition2_frame <- data.frame(condition2)
condition2_sorted <- condition2_frame[ , order(names(condition2_frame))]
# replace zeros in dataframe
condition2_sorted %>% remove_rownames %>% column_to_rownames(var="a_0gene") -> condition2_final
condition1_sorted %>% remove_rownames %>% column_to_rownames(var="a_0gene") -> condition1_final
# rename genes
# labeled_data <- rename_genes(data = data, id = "ensembl", change = "names")

(condition <- factor(c(rep("ctl", 4), rep("exp", 80))))
(samples <- factor(substr(colnames(condition1), 1, 1)))
# Convert to matrix
countdata1 <- as.matrix(condition1_final)
countdata2 <- as.matrix(condition2_final)
countdata1[is.na(countdata1)] <- 0
countdata2[is.na(countdata2)] <- 0
# make DESeq2 dataframe
(coldata1 <- data.frame(row.names=colnames(countdata1), condition))
(coldata2 <- data.frame(row.names=colnames(countdata2), condition))
dds1 <- DESeqDataSetFromMatrix(countData=countdata1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=countdata2, colData=coldata2, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

# plot dispersions
# pdf("qc-dispersions_condition1.pdf", pointsize=.5)
# plotDispEsts(dds1, main="Dispersion plot")
# dev.off()
# 
# pdf("qc-dispersions_condition2.pdf", pointsize=.5)
# plotDispEsts(dds2, main="Dispersion plot") 
# dev.off()

# Regularized log transformation for clustering/heatmaps, etc
# rld1 <- rlogTransformation(dds1)
# rld2 <- rlogTransformation(dds2)
# saveRDS(rld1, "./rlogtransform1.rds")
# saveRDS(rld2, "./rlogtransform2.rds")

rld1 <- readRDS("data/rlogtransform1.rds")
rld2 <- readRDS("data/rlogtransform2.rds")

hist(assay(rld1))
hist(assay(rld2))


# (mycols <- c("#FFD700", "#4682B4"))
# reorder


############
# DISTANCE HEATMAP
############


# Sample distance heatmap
sampleDists1 <- as.matrix(dist(t(assay(rld1))))
colnames(sampleDists1) <- substr(colnames(sampleDists1), 1, 1)
rownames(sampleDists1) <- substr(rownames(sampleDists1), 1, 1)

# color dendogram branches
hcluster = hclust(dist(t(assay(rld1))), method ="ward.D")
dend1 <- as.dendrogram(hcluster)
library(dendextend)
(mycols <- sample(heat.colors(21)[1:length(unique(colnames(sampleDists1)))]))
dend1 <- color_branches(dend1, k = 21, col = mycols)
col_labels <- get_leaves_branches_col(dend1)
# But due to the way heatmap.2 works - we need to fix it to be in the 
# order of the data!    
col_labels <- col_labels[order(order.dendrogram(dend1))]

pdf("qc-heatmap-samples1.pdf", pointsize=.5)
heatmap.2(as.matrix(sampleDists1), key=F, trace="none",
          col=colorpanel(100, "#FFD700", "#4682B4"),
          cluster_columns = FALSE,
          # dendrogram="row",
          Rowv = dend1,
          # Colv = dend1,
          ColSideColors=mycols[samples],
          RowSideColors=col_labels,
          margin=c(10, 10), main="Sample Distance Matrix Condition 1")
dev.off()


sampleDists2 <- as.matrix(dist(t(assay(rld2))))
colnames(sampleDists2) <- substr(colnames(sampleDists2), 1, 1)
rownames(sampleDists2) <- substr(rownames(sampleDists2), 1, 1)
(mycols <- sample(rainbow(21)[1:length(unique(colnames(sampleDists2)))]))
pdf("qc-heatmap-samples2.pdf", pointsize=.5)
heatmap.2(as.matrix(sampleDists2), key=F, trace="none",
          col=colorpanel(100, "#FFD700", "#4682B4"),
          cluster_columns = FALSE,
          ColSideColors=mycols[samples], RowSideColors=mycols[samples],
          margin=c(10, 10), main="Sample Distance Matrix Condition 2")
dev.off()

# pdf("qc-pca1.pdf", pointsize=10)
# rld_pca(rld1, intgroup="condition", xlim=c(-75, 35))
# dev.off()
######################################################################################

##########
# PCA of samples
###########

pdf("qc-pca_labs.pdf", pointsize=10)
pca1 <- rld_pca(rld1, intgroup="condition", xlim=c(-20, 20))
dev.off()

pdf("qc-pca-scaled.pdf", pointsize=10)
rld_pca(rld1, intgroup="condition", xlim=c(-15, 15))
dev.off()

pdf("qc-pca-scaled_nolabs.pdf", pointsize=10)
rld_pca(rld1, intgroup="condition", xlim=c(-15, 15))
dev.off()

pdf("qc-pca-scaled_nolabs2.pdf", pointsize=10)
pca2 <- rld_pca(rld2, intgroup="condition", xlim=c(-15, 15))
dev.off()
########
# silhouette analysis
########

# silhouette for condition 1
silhouette1 <- silhouette(pca1)

# silhouette for condition 2
rownames(pca2$x) <- substr(rownames(pca2$x), 4, 6)
rownames(pca2$x)[1] <- "A_1"
rownames(pca2$x)[2] <- "A_2"
rownames(pca2$x)[3] <- "A_3"
rownames(pca2$x)[4] <- "A_4"

silhouette2 <- silhouette(pca2)


