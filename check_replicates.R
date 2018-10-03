#####
# EXECUTE
#####
library(dplyr)
library(tibble)
library(RColorBrewer)
library(gplots)
library(genefilter)
library(calibrate)
library(DESeq2)
# source original function list 
source("./helper_functions.R")

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
# Convert to matrix
countdata1 <- as.matrix(condition1_averaged)
countdata2 <- as.matrix(condition2_averaged)

# make DESeq2 dataframe
(coldata1 <- data.frame(row.names=colnames(countdata1), condition))
(coldata2 <- data.frame(row.names=colnames(countdata2), condition))
dds1 <- DESeqDataSetFromMatrix(countData=countdata1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=countdata2, colData=coldata2, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

# plot dispersions
png("qc-dispersions_condition1.png", 1000, 1000, pointsize=20)
plotDispEsts(dds1, main="Dispersion plot")
dev.off()

png("qc-dispersions_condition2.png", 1000, 1000, pointsize=20)
plotDispEsts(dds2, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld1 <- rlogTransformation(dds1)
rld2 <- rlogTransformation(dds2)
head(assay(rld))
hist(assay(rld))

(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists1 <- as.matrix(dist(t(assay(rld1))))
png("qc-heatmap-samples1.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists1), key=F, trace="none",
          col=colorpanel(100, "red", "green"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()


rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca1.png", 1000, 1000, pointsize=20)
rld_pca(rld1, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()


# get differential expression results
res <- results(dds2)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)