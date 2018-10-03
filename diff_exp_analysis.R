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

# rename columns to look nice
colnames(countdata) <- gsub("\\.fastq.trim.sam.counts", "", colnames(countdata))
colnames(countdata) <- gsub("X", "", colnames(countdata))
colnames(a_data) <- gsub("\\.fastq.trim.sam.counts", "", colnames(a_data))
colnames(a_data) <- gsub("X", "", colnames(a_data))

# merege A and allcounts
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

# rename rownames from numbers to gene names (for R object)
condition2_sorted %>% remove_rownames %>% column_to_rownames(var="a_0gene") -> condition2_final
condition1_sorted %>% remove_rownames %>% column_to_rownames(var="a_0gene") -> condition1_final
# rename genes
# labeled_data <- rename_genes(data = data, id = "ensembl", change = "names")

# Remove NAs from both data frames
# replace zeros in dataframe
condition1_final[is.na(condition1_final)] <- 0
condition2_final[is.na(condition2_final)] <- 0

# average across columns
condition1_averaged <- t(apply(condition1_final, 1, tapply, gl(21, 4), mean))

# average across columns
condition2_averaged <- t(apply(condition2_final, 1, tapply, gl(21, 4), mean))

# rename colnames 
alpha <- LETTERS[1:21]
for (i in 1:21){
  colnames(condition1_averaged)[i] <- alpha[i]
}

for (i in 1:21){
  colnames(condition2_averaged)[i] <- paste("2", alpha[i], sep="")
}



#########################################################

############
# RUN DESEQ2
############

# set up experimental parameters
(condition <- factor(c(rep("ctl", 1), rep("exp", 20))))

# Convert to matrix
countdata1 <- as.matrix(condition1_averaged)
countdata1 <- apply(countdata1, 2, as.integer)
rownames(countdata1) <- c(rownames(condition1_averaged))

countdata2 <- as.matrix(condition2_averaged)
countdata2 <- apply(countdata2, 2, as.integer)
rownames(countdata2) <- c(rownames(condition2_averaged))

# non-normalized
png("non_normalized1.png", 2000,2000, pointsize=20)
heatmap(head(countdata1, 10000), Colv=NA)
dev.off()

png("non_normalized2.png", 2000,2000, pointsize=20)
heatmap(head(countdata2, 10000), Colv=NA)
dev.off()

# make DESeq2 dataframe
(coldata1 <- data.frame(row.names=colnames(countdata1), condition))
(coldata2 <- data.frame(row.names=colnames(countdata2), condition))
dds1 <- DESeqDataSetFromMatrix(countData=countdata1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=countdata2, colData=coldata2, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)


# get differential expression results
res1 <- results(dds1)
table(res1$pvalue<0.01)
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

res2 <- results(dds2)
table(res2$padj<0.05)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

# make heatmaps
hmap1 <- as.matrix(resdata1[1:10000, -(1:7)])
rownames(hmap1) <- resdata1[1:10000 ,1]
png("heatmap_condition1_raw.png", 2000,2000, pointsize=20)
heatmap(hmap1, Colv = NA)
dev.off()

hmap2 <- as.matrix(resdata2[1:100, -(1:7)])
rownames(hmap2) <- resdata2[1:100 ,1]
png("heatmap_condition2_raw.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()

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


######
# Generate exploratory PCA plot
######

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

###############
# CONDITION 1
###############

# subset data


countdata1 <- as.data.frame(countdata1)

#################################################
countdata1 %>% dplyr::select(B, G, L, Q) -> fifteen
rownames(fifteen) <- rownames(countdata1)
fifteen <- apply(fifteen, 2, as.numeric)

png("fifteen_minutes.png", 2000, 2000, pointsize=20)
heatmap(head(fifteen, 10000), Colv=NA)
dev.off()
################################################
countdata1 %>% dplyr::select(C, H, M, R) -> thirty
thirty <- apply(thirty, 2, as.numeric)

png("thirty_minutes.png", 2000, 2000, pointsize=20)
heatmap(head(thirty, 10000), Colv=NA)
dev.off()

################################################
countdata1 %>% dplyr::select(D, I, N, S) -> hour1
hour1 <- apply(hour1, 2, as.numeric)

png("hour1.png", 2000, 2000, pointsize=20)
heatmap(head(hour1, 10000), Colv=NA)
dev.off()

##################################################
countdata1 %>% dplyr::select(E, J, O, T) -> hour2
hour2 <- apply(hour2, 2, as.numeric)

png("hour2.png", 2000, 2000, pointsize=20)
heatmap(head(hour2, 10000), Colv=NA)
dev.off()

#################################################
countdata1 %>% dplyr::select(F, K, P, U) -> hour6
hour6 <- apply(hour6, 2, as.numeric)

png("hour6.png", 2000, 2000, pointsize=20)
heatmap(head(hour6, 10000), Colv=NA)
dev.off()

###############
# CONDITION 2
###############
countdata2 <- as.data.frame(countdata2)

#################################################
countdata2 %>% dplyr::select("2B", "2G", "2L", "2Q") -> fifteen2
fifteen2 <- apply(fifteen2, 2, as.numeric)

png("2_fifteen_minutes.png", 2000, 2000, pointsize=20)
heatmap(head(fifteen2, 10000), Colv=NA)
dev.off()
################################################
countdata2 %>% dplyr::select("2C", "2H", "2M", "2R") -> thirty2
thirty2 <- apply(thirty2, 2, as.numeric)

png("2_thirty_minutes.png", 2000, 2000, pointsize=20)
heatmap(head(thirty2, 10000), Colv=NA)
dev.off()

################################################
countdata2 %>% dplyr::select("2D", "2I", "2N", "2S") -> hour1_2
hour1_2 <- apply(hour1_2, 2, as.numeric)

png("2_hour1.png", 2000, 2000, pointsize=20)
heatmap(head(hour1_2, 10000), Colv=NA)
dev.off()

##################################################
countdata2 %>% dplyr::select("2E", "2J", "2O", "2T") -> hour2_2
hour2_2 <- apply(hour2_2, 2, as.numeric)

png("2_hour2.png", 2000, 2000, pointsize=20)
heatmap(head(hour2_2, 10000), Colv=NA)
dev.off()

#################################################
countdata2 %>% dplyr::select("2F", "2K", "2P", "2U") -> hour6_2
hour6_2 <- apply(hour6_2, 2, as.numeric)

png("2_hour6.png", 2000, 2000, pointsize=20)
heatmap(head(hour6_2, 10000), Colv=NA)
dev.off()

#############
# ALL SAMPLES HEAT MAPS
############
png("non_normalized1.png", 2000,2000, pointsize=20)
heatmap(head(countdata1, 10000), Colv=NA)
dev.off()

png("non_normalized2.png", 2000,2000, pointsize=20)
heatmap(head(countdata2, 10000), Colv=NA)
dev.off()


#######
# Generate differential expression analysis of subset 1
#######
######################################
# subset data 
(condition <- factor(c(rep("ctl", 1), rep("exp", 3))))
# make DESeq2 dataframe
(coldata1 <- data.frame(row.names=colnames(fifteen), condition))
(coldata2 <- data.frame(row.names=colnames(fifteen2), condition))
dds1 <- DESeqDataSetFromMatrix(countData=fifteen, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=fifteen2, colData=coldata2, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

# get differential expression results
res1 <- results(dds1)
table(res1$padj<0.05)
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)


res2 <- results(dds2)
table(res2$padj<0.05)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

##########################





#########
# 15 min
#########
fifteen <- as.data.frame(fifteen)

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe
(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))


countdata %>% dplyr::select(B_1, B_2, B_3, B_4, G_1, G_2, G_3, G_4) -> pairwise1
countdata %>% dplyr::select(B_1, B_2, B_3, B_4, L_1, L_2, L_3, L_4) -> pairwise2
countdata %>% dplyr::select(B_1, B_2, B_3, B_4, Q_1, Q_2, Q_3, Q_4) -> pairwise3

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
a <- table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("heatmap_condition2_raw.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
png("./BvsL_MA.png", 2000,2000, pointsize=60)
plotMA(res2, ylim=c(-2,2))
dev.off()
table(res2$padj<0.01)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
table(res3$padj<0.01)
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)


#######
# 30 min
#######
thirty <- as.data.frame(thirty)

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select(C_1, C_2, C_3, C_4, H_1, H_2, H_3, H_4) -> pairwise1
countdata %>% dplyr::select(C_1, C_2, C_3, C_4, M_1, M_2, M_3, M_4) -> pairwise2
countdata %>% dplyr::select(C_1, C_2, C_3, C_4, R_1, R_2, R_3, R_4) -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))



dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
a <- table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("condition1_thirty_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
png("./CvsM_MA.png", 2000,2000, pointsize=60)
plotMA(res2, ylim=c(-2,2))
dev.off()
table(res2$padj<0.01)
a <- table(res2$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)
resdata2_final <- head(resdata2, number_to_keep)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
a <- table(res3$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"

resdata3_final <- head(resdata3, number_to_keep)

res_data1_final <- as.data.frame(resdata1_final[, -(2:7)])
res_data2_final <- as.data.frame(resdata2_final[, -(2:11)])
res_data3_final <- as.data.frame(resdata3_final[, -(2:11)])


middle <- merge(res_data1_final, res_data2_final, by="Gene", no.dups=TRUE)
final <- merge(middle, res_data3_final, by="Gene", no.dups=TRUE)
final <- na.omit(final)
final <- final[order(final$padj), ]
hmap2 <- as.matrix(final[, -1])
rownames(hmap2) <- final[ ,1]
png("thirty_min_normalized_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()


#######
# 60 MIN
#######
hour1 <- as.data.frame(hour1)

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select(D_1, D_2, D_3, D_4, I_1, I_2, I_3, I_4) -> pairwise1
countdata %>% dplyr::select(D_1, D_2, D_3, D_4, N_1, N_2, N_3, N_4) -> pairwise2
countdata %>% dplyr::select(D_1, D_2, D_3, D_4, S_1, S_2, S_3, S_4) -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))



dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
a <- table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("condition1_thirty_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
png("./CvsM_MA.png", 2000,2000, pointsize=60)
plotMA(res2, ylim=c(-2,2))
dev.off()
table(res2$padj<0.01)
a <- table(res2$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)
resdata2_final <- head(resdata2, number_to_keep)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
a <- table(res3$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"

resdata3_final <- head(resdata3, number_to_keep)

res_data1_final <- as.data.frame(resdata1_final[, -(2:7)])
res_data2_final <- as.data.frame(resdata2_final[, -(2:11)])
res_data3_final <- as.data.frame(resdata3_final[, -(2:11)])


middle <- merge(res_data1_final, res_data2_final, by="Gene", no.dups=TRUE)
final <- merge(middle, res_data3_final, by="Gene", no.dups=TRUE)
final <- na.omit(final)
final <- final[order(final$padj), ]
hmap2 <- as.matrix(final[, -1])
rownames(hmap2) <- final[ ,1]
png("thirty_min_normalized_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()

######
# 120 min
######

hour2 <- as.data.frame(hour2)

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select(E_1, E_2, E_3, E_4, J_1, J_2, J_3, J_4) -> pairwise1
countdata %>% dplyr::select(E_1, E_2, E_3, E_4, O_1, O_2, O_3, O_4) -> pairwise2
countdata %>% dplyr::select(E_1, E_2, E_3, E_4, T_1, T_2, T_3, T_4) -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))



dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
a <- table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("condition1_thirty_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
png("./CvsM_MA.png", 2000,2000, pointsize=60)
plotMA(res2, ylim=c(-2,2))
dev.off()
table(res2$padj<0.01)
a <- table(res2$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)
resdata2_final <- head(resdata2, number_to_keep)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
a <- table(res3$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"

resdata3_final <- head(resdata3, number_to_keep)

res_data1_final <- as.data.frame(resdata1_final[, -(2:7)])
res_data2_final <- as.data.frame(resdata2_final[, -(2:11)])
res_data3_final <- as.data.frame(resdata3_final[, -(2:11)])


middle <- merge(res_data1_final, res_data2_final, by="Gene", no.dups=TRUE)
final <- merge(middle, res_data3_final, by="Gene", no.dups=TRUE)
final <- na.omit(final)
final <- final[order(final$padj), ]
hmap2 <- as.matrix(final[, -1])
rownames(hmap2) <- final[ ,1]
png("thirty_min_normalized_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()


#######
# 6 hours
#######

hour6 <- as.data.frame(hour6)

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select(F_1, F_2, F_3, F_4, K_1, K_2, K_3, K_4) -> pairwise1
countdata %>% dplyr::select(F_1, F_2, F_3, F_4, P_1, P_2, P_3, P_4) -> pairwise2
countdata %>% dplyr::select(F_1, F_2, F_3, F_4, U_1, U_2, U_3, U_4) -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))



dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
a <- table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("condition1_thirty_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
png("./CvsM_MA.png", 2000,2000, pointsize=60)
plotMA(res2, ylim=c(-2,2))
dev.off()
table(res2$padj<0.01)
a <- table(res2$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)
resdata2_final <- head(resdata2, number_to_keep)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
a <- table(res3$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"

resdata3_final <- head(resdata3, number_to_keep)

res_data1_final <- as.data.frame(resdata1_final[, -(2:7)])
res_data2_final <- as.data.frame(resdata2_final[, -(2:11)])
res_data3_final <- as.data.frame(resdata3_final[, -(2:11)])


middle <- merge(res_data1_final, res_data2_final, by="Gene", no.dups=TRUE)
final <- merge(middle, res_data3_final, by="Gene", no.dups=TRUE)
final <- na.omit(final)
final <- final[order(final$padj), ]
hmap2 <- as.matrix(final[, -1])
rownames(hmap2) <- final[ ,1]
png("thirty_min_normalized_heatmap.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()
# visualize differential expression results 

############################################### CONDITION 2 ##############################################

#########
# 15 MIN
#########

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_B_1", "2_B_2", "2_B_3", "2_B_4", "2_G_1", "2_G_2", "2_G_3", "2_G_4") -> pairwise1
countdata %>% dplyr::select("2_B_1", "2_B_2", "2_B_3", "2_B_4", "2_L_1", "2_L_2", "2_L_3", "2_L_4") -> pairwise2
countdata %>% dplyr::select("2_B_1", "2_B_2", "2_B_3", "2_B_4", "2_Q_1", "2_Q_2", "2_Q_3", "2_Q_4") -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("heatmap_condition2_raw.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
plotMA(res2, ylim=c(-2,2))
table(res2$padj<0.01)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
table(res3$padj<0.01)
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)


######
# 30 MIN
######

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
modified <- factor(c(rep("ctl", 4), rep("exp", 3)))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_C_1", "2_C_2", "2_C_3", "2_C_4", "2_H_1", "2_H_2", "2_H_3", "2_H_4") -> pairwise1
countdata %>% dplyr::select("2_C_1", "2_C_2", "2_C_3", "2_C_4", "2_M_2", "2_M_3", "2_M_4") -> pairwise2
countdata %>% dplyr::select("2_C_1", "2_C_2", "2_C_3", "2_C_4", "2_R_1", "2_R_2", "2_R_3", "2_R_4") -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), modified))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~modified)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("heatmap_condition2_raw.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
plotMA(res2, ylim=c(-2,2))
table(res2$padj<0.01)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
table(res3$padj<0.01)
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

#####
# 60 min
####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_D_1", "2_D_2", "2_D_3", "2_D_4", "2_I_1", "2_I_2", "2_I_3", "2_I_4") -> pairwise1
countdata %>% dplyr::select("2_D_1", "2_D_2", "2_D_3", "2_D_4", "2_N_1", "2_N_2", "2_N_3", "2_N_4") -> pairwise2
countdata %>% dplyr::select("2_D_1", "2_D_2", "2_D_3", "2_D_4", "2_S_1", "2_S_2", "2_S_3", "2_S_4") -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("heatmap_condition2_raw.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
plotMA(res2, ylim=c(-2,2))
table(res2$padj<0.01)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
table(res3$padj<0.01)
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

####
# 120 MIN
####
(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_E_1", "2_E_2", "2_E_3", "2_E_4", "2_J_1", "2_J_2", "2_J_3", "2_J_4") -> pairwise1
countdata %>% dplyr::select("2_E_1", "2_E_2", "2_E_3", "2_E_4", "2_O_1", "2_O_2", "2_O_3", "2_O_4") -> pairwise2
countdata %>% dplyr::select("2_E_1", "2_E_2", "2_E_3", "2_E_4", "2_T_1", "2_T_2", "2_T_3", "2_T_4") -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("heatmap_condition2_raw.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
plotMA(res2, ylim=c(-2,2))
table(res2$padj<0.01)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
table(res3$padj<0.01)
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

#####
# 360 min
#####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_F_1", "2_F_2", "2_F_3", "2_F_4", "2_K_1", "2_K_2", "2_K_3", "2_K_4") -> pairwise1
countdata %>% dplyr::select("2_F_1", "2_F_2", "2_F_3", "2_F_4", "2_P_1", "2_P_2", "2_P_3", "2_P_4") -> pairwise2
countdata %>% dplyr::select("2_F_1", "2_F_2", "2_F_3", "2_F_4", "2_U_1", "2_U_2", "2_U_3", "2_U_4") -> pairwise3

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)


res1 <- results(dds1)
plotMA(res1, ylim=c(-2,2))
table(res1$padj<0.01)
number_to_keep <- a[2]
## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)
resdata1_final <- head(resdata1, number_to_keep)

hmap2 <- as.matrix(resdata1_final[, -(1:7)])
rownames(hmap2) <- resdata1_final[ ,1]
png("heatmap_condition2_raw.png", 2000,2000, pointsize=20)
heatmap(hmap2, Colv = NA)
dev.off()



res2 <- results(dds2)
plotMA(res2, ylim=c(-2,2))
table(res2$padj<0.01)
## Order by adjusted p-value
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)

res3 <- results(dds3)
plotMA(res3, ylim=c(-2,2))
table(res3$padj<0.01)
## Order by adjusted p-value
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)


########################################################## Time analysis ##########################################################################
########################################################## Condition 2   ##########################################################################
#####
# 25 micromolar
#####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_B_1", "2_B_2", "2_B_3", "2_B_4", "2_G_1", "2_G_2", "2_G_3", "2_G_4") -> pairwise1
countdata %>% dplyr::select("2_C_1", "2_C_2", "2_C_3", "2_C_4", "2_H_1", "2_H_2", "2_H_3", "2_H_4") -> pairwise2
countdata %>% dplyr::select("2_D_1", "2_D_2", "2_D_3", "2_D_4", "2_I_1", "2_I_2", "2_I_3", "2_I_4") -> pairwise3
countdata %>% dplyr::select("2_E_1", "2_E_2", "2_E_3", "2_E_4", "2_J_1", "2_J_2", "2_J_3", "2_J_4") -> pairwise4
countdata %>% dplyr::select("2_F_1", "2_F_2", "2_F_3", "2_F_4", "2_K_1", "2_K_2", "2_K_3", "2_K_4") -> pairwise5

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))
(coldata4 <- data.frame(row.names=colnames(pairwise4), condition))
(coldata5 <- data.frame(row.names=colnames(pairwise5), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds4 <- DESeqDataSetFromMatrix(countData=pairwise4, colData=coldata4, design=~condition)
dds5 <- DESeqDataSetFromMatrix(countData=pairwise5, colData=coldata5, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)
dds4 <- DESeq(dds4)
dds5 <- DESeq(dds5)

res1 <- results(dds1)
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

res2 <- results(dds2)
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)


res3 <- results(dds3)
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

res4 <- results(dds4)
res4 <- res4[order(res4$padj), ]
## Merge with normalized count data
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata4)[1] <- "Gene"
head(resdata4)

res5 <- results(dds5)
res5 <- res5[order(res5$padj), ]
## Merge with normalized count data
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata5)[1] <- "Gene"
head(resdata5)



#####
# 50 micromolar
#####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_B_1", "2_B_2", "2_B_3", "2_B_4", "2_L_1", "2_L_2", "2_L_3", "2_L_4") -> pairwise1
countdata %>% dplyr::select("2_C_1", "2_C_2", "2_C_3", "2_C_4", "2_M_1", "2_M_2", "2_M_3", "2_M_4") -> pairwise2
countdata %>% dplyr::select("2_D_1", "2_D_2", "2_D_3", "2_D_4", "2_N_1", "2_N_2", "2_N_3", "2_N_4") -> pairwise3
countdata %>% dplyr::select("2_E_1", "2_E_2", "2_E_3", "2_E_4", "2_O_1", "2_O_2", "2_O_3", "2_O_4") -> pairwise4
countdata %>% dplyr::select("2_F_1", "2_F_2", "2_F_3", "2_F_4", "2_P_1", "2_P_2", "2_P_3", "2_P_4") -> pairwise5

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))
(coldata4 <- data.frame(row.names=colnames(pairwise4), condition))
(coldata5 <- data.frame(row.names=colnames(pairwise5), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds4 <- DESeqDataSetFromMatrix(countData=pairwise4, colData=coldata4, design=~condition)
dds5 <- DESeqDataSetFromMatrix(countData=pairwise5, colData=coldata5, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)
dds4 <- DESeq(dds4)
dds5 <- DESeq(dds5)

res1 <- results(dds1)
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

res2 <- results(dds2)
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)


res3 <- results(dds3)
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

res4 <- results(dds4)
res4 <- res4[order(res4$padj), ]
## Merge with normalized count data
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata4)[1] <- "Gene"
head(resdata4)

res5 <- results(dds5)
res5 <- res5[order(res5$padj), ]
## Merge with normalized count data
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata5)[1] <- "Gene"
head(resdata5)

####
# 100 micromolar
####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("2_B_1", "2_B_2", "2_B_3", "2_B_4", "2_Q_1", "2_Q_2", "2_Q_3", "2_Q_4") -> pairwise1
countdata %>% dplyr::select("2_C_1", "2_C_2", "2_C_3", "2_C_4", "2_R_1", "2_R_2", "2_R_3", "2_R_4") -> pairwise2
countdata %>% dplyr::select("2_D_1", "2_D_2", "2_D_3", "2_D_4", "2_S_1", "2_S_2", "2_S_3", "2_S_4") -> pairwise3
countdata %>% dplyr::select("2_E_1", "2_E_2", "2_E_3", "2_E_4", "2_T_1", "2_T_2", "2_T_3", "2_T_4") -> pairwise4
countdata %>% dplyr::select("2_F_1", "2_F_2", "2_F_3", "2_F_4", "2_U_1", "2_U_2", "2_U_3", "2_U_4") -> pairwise5

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))
(coldata4 <- data.frame(row.names=colnames(pairwise4), condition))
(coldata5 <- data.frame(row.names=colnames(pairwise5), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds4 <- DESeqDataSetFromMatrix(countData=pairwise4, colData=coldata4, design=~condition)
dds5 <- DESeqDataSetFromMatrix(countData=pairwise5, colData=coldata5, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)
dds4 <- DESeq(dds4)
dds5 <- DESeq(dds5)

res1 <- results(dds1)
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

res2 <- results(dds2)
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)


res3 <- results(dds3)
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

res4 <- results(dds4)
res4 <- res4[order(res4$padj), ]
## Merge with normalized count data
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata4)[1] <- "Gene"
head(resdata4)

res5 <- results(dds5)
res5 <- res5[order(res5$padj), ]
## Merge with normalized count data
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata5)[1] <- "Gene"
head(resdata5)


########################################################### CONDITION 1 #############################################################


#####
# 25 micromolar
#####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("B_1", "B_2", "B_3", "B_4", "G_1", "G_2", "G_3", "G_4") -> pairwise1
countdata %>% dplyr::select("C_1", "C_2", "C_3", "C_4", "H_1", "H_2", "H_3", "H_4") -> pairwise2
countdata %>% dplyr::select("D_1", "D_2", "D_3", "D_4", "I_1", "I_2", "I_3", "I_4") -> pairwise3
countdata %>% dplyr::select("E_1", "E_2", "E_3", "E_4", "J_1", "J_2", "J_3", "J_4") -> pairwise4
countdata %>% dplyr::select("F_1", "F_2", "F_3", "F_4", "K_1", "K_2", "K_3", "K_4") -> pairwise5

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))
(coldata4 <- data.frame(row.names=colnames(pairwise4), condition))
(coldata5 <- data.frame(row.names=colnames(pairwise5), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds4 <- DESeqDataSetFromMatrix(countData=pairwise4, colData=coldata4, design=~condition)
dds5 <- DESeqDataSetFromMatrix(countData=pairwise5, colData=coldata5, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)
dds4 <- DESeq(dds4)
dds5 <- DESeq(dds5)

res1 <- results(dds1)
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

res2 <- results(dds2)
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)


res3 <- results(dds3)
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

res4 <- results(dds4)
res4 <- res4[order(res4$padj), ]
## Merge with normalized count data
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata4)[1] <- "Gene"
head(resdata4)

res5 <- results(dds5)
res5 <- res5[order(res5$padj), ]
## Merge with normalized count data
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata5)[1] <- "Gene"
head(resdata5)



#####
# 50 micromolar
#####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("B_1", "B_2", "B_3", "B_4", "L_1", "L_2", "L_3", "L_4") -> pairwise1
countdata %>% dplyr::select("C_1", "C_2", "C_3", "C_4", "M_1", "M_2", "M_3", "M_4") -> pairwise2
countdata %>% dplyr::select("D_1", "D_2", "D_3", "D_4", "N_1", "N_2", "N_3", "N_4") -> pairwise3
countdata %>% dplyr::select("E_1", "E_2", "E_3", "E_4", "O_1", "O_2", "O_3", "O_4") -> pairwise4
countdata %>% dplyr::select("F_1", "F_2", "F_3", "F_4", "P_1", "P_2", "P_3", "P_4") -> pairwise5

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))
(coldata4 <- data.frame(row.names=colnames(pairwise4), condition))
(coldata5 <- data.frame(row.names=colnames(pairwise5), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds4 <- DESeqDataSetFromMatrix(countData=pairwise4, colData=coldata4, design=~condition)
dds5 <- DESeqDataSetFromMatrix(countData=pairwise5, colData=coldata5, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)
dds4 <- DESeq(dds4)
dds5 <- DESeq(dds5)

res1 <- results(dds1)
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

res2 <- results(dds2)
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)


res3 <- results(dds3)
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

res4 <- results(dds4)
res4 <- res4[order(res4$padj), ]
## Merge with normalized count data
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata4)[1] <- "Gene"
head(resdata4)

res5 <- results(dds5)
res5 <- res5[order(res5$padj), ]
## Merge with normalized count data
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata5)[1] <- "Gene"
head(resdata5)

####
# 100 micromolar
####

(condition <- factor(c(rep("ctl", 4), rep("exp", 4))))
# make DESeq2 dataframe

countdata %>% dplyr::select("B_1", "B_2", "B_3", "B_4", "Q_1", "Q_2", "Q_3", "Q_4") -> pairwise1
countdata %>% dplyr::select("C_1", "C_2", "C_3", "C_4", "R_1", "R_2", "R_3", "R_4") -> pairwise2
countdata %>% dplyr::select("D_1", "D_2", "D_3", "D_4", "S_1", "S_2", "S_3", "S_4") -> pairwise3
countdata %>% dplyr::select("E_1", "E_2", "E_3", "E_4", "T_1", "T_2", "T_3", "T_4") -> pairwise4
countdata %>% dplyr::select("F_1", "F_2", "F_3", "F_4", "U_1", "U_2", "U_3", "U_4") -> pairwise5

(coldata1 <- data.frame(row.names=colnames(pairwise1), condition))
(coldata2 <- data.frame(row.names=colnames(pairwise2), condition))
(coldata3 <- data.frame(row.names=colnames(pairwise3), condition))
(coldata4 <- data.frame(row.names=colnames(pairwise4), condition))
(coldata5 <- data.frame(row.names=colnames(pairwise5), condition))

dds1 <- DESeqDataSetFromMatrix(countData=pairwise1, colData=coldata1, design=~condition)
dds2 <- DESeqDataSetFromMatrix(countData=pairwise2, colData=coldata2, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=pairwise3, colData=coldata3, design=~condition)
dds4 <- DESeqDataSetFromMatrix(countData=pairwise4, colData=coldata4, design=~condition)
dds5 <- DESeqDataSetFromMatrix(countData=pairwise5, colData=coldata5, design=~condition)
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)
dds4 <- DESeq(dds4)
dds5 <- DESeq(dds5)

res1 <- results(dds1)
res1 <- res1[order(res1$padj), ]
## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

res2 <- results(dds2)
res2 <- res2[order(res2$padj), ]
## Merge with normalized count data
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata2)[1] <- "Gene"
head(resdata2)


res3 <- results(dds3)
res3 <- res3[order(res3$padj), ]
## Merge with normalized count data
resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata3)[1] <- "Gene"
head(resdata3)

res4 <- results(dds4)
res4 <- res4[order(res4$padj), ]
## Merge with normalized count data
resdata4 <- merge(as.data.frame(res4), as.data.frame(counts(dds4, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata4)[1] <- "Gene"
head(resdata4)

res5 <- results(dds5)
res5 <- res5[order(res5$padj), ]
## Merge with normalized count data
resdata5 <- merge(as.data.frame(res5), as.data.frame(counts(dds5, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata5)[1] <- "Gene"
head(resdata5)





























#######
# Generate differential expression analysis of subset 1
#######

# This analysis is inhibitor concentration



#######
# Generate differential expression analysis of subset 1
#######




#######
# Determine peak expression at no inhibitor condition
#######
countdata1_ugh <- apply(countdata1, 2, as.numeric)
colMeans(countdata1_ugh)


#######
# Make line graphs of differentially expressed genes
#######

# get ids for 15 minutes 
diff_ids <- data.frame(gene = head(resdata1$Gene, 5), sig = head(resdata1$padj, 5))
fifteen_id <- c("B", "G", "L", "Q")
# countdata1["gene"] <- rownames(countdata1)
diff_frame_15 <- merge(diff_ids, countdata1, by.x = "gene")
# diff_frame_15 <- subset(countdata1, rownames(condition1_final) %in% diff_ids)
diff_frame_15 %>% dplyr::select(gene,fifteen_id) -> diff_frame_15

# diff_frame_15t <- data.frame(t(diff_frame_15))
# ["gene"] <- rownames(diff_frame_15)
diff_frame_15.molten <- melt(diff_frame_15, id.vars="gene", value.name="reads", variable.name="inhibitor")
ggplot(diff_frame_15.molten, aes(x= inhibitor, y=reads, color=factor(gene))) + geom_path(aes_string(group="gene")) + geom_point() + fte_theme()


######
# Make line graphs of single genes with respect to time
######

allcounts %>% filter(rn == "ENSMUSG00000026069.15") -> il1
View(il1)
il1.molten <- melt(il1, id.vars = "rn", value.name = "reads")
order(il1.molten)
