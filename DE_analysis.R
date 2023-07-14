# Download necessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")

## Transcript level analysis

# Load sample file information
setwd("~/RNA_122022")
csvfile <- file.path("SampleSSvsPP.csv")
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)
coldata
str(coldata)

# Convert Condition from character to factor
coldata$Condition <- as.factor(coldata$Condition)
str(coldata)

# Quantification file based on SampleTable information
files <- file.path("quantSSvsPP", coldata$Sample, "quantFinalcor.sf")
files

# Assign names to files
names(files) <- coldata$Sample
library("readr")

# Import the necessary quantification data for DESeq2 using the tximport function
library("tximport")
library("jsonlite")
txi.tr <- tximport(files, type = "salmon", txOut = TRUE)
txi.tr
?tximport

# Check count per sample
head(txi.tr$counts)

# Load DESeq2 library
library("DESeq2")

# Create DESeqDataSet object
ddsTxi.tr <- DESeqDataSetFromTximport(txi.tr, colData = coldata, design = ~ Condition)

#### Count data (independent of expression analysis) for k-means clustering
# No filter on read count

ddsTxi.tr <- DESeq(ddsTxi.tr)
ddsTxi.tr

# Retrieve the normalized count matrix
# Normalized count
mat <- counts(ddsTxi.tr, normalized = TRUE)
mat
write.csv(mat, file = "norm_countsPPvsSS.csv")

# Unnormalized count
mat0 <- counts(ddsTxi.tr, normalized = FALSE)
mat0
write.csv(mat0, file = "unnorm_countsPPvsSS.csv")


#################################################################

## Set reference level

ddsTxi.tr$Condition <- relevel(ddsTxi.tr$Condition, "C3")

levels(ddsTxi.tr$Condition)
colData(ddsTxi.tr)

## Perform differential expression analysis

dds.tr <- DESeq(ddsTxi.tr)
res.tr <- results(dds.tr)
res.tr
summary(res.tr)
resultsNames(dds.tr)

# Extract the general results and order by p-value
resC3vsC4allOrdered <- res.tr[order(res.tr$pvalue), ]
resC3vsC4allOrdered

# Write the general results
write.csv(as.data.frame(resC3vsC4allOrdered), file = "resC3vsC4ALLOrdered.csv")

## Extract specific results

resC3_vs_C4 <- results(dds.tr, name = "Condition_C4_vs_C3", lfcThreshold = 1, alpha = 0.05)
resC3_vs_C4
summary(resC3_vs_C4)

# Order results by p-value
resC3_vs_C4Ordered <- resC3_vs_C4[order(resC3_vs_C4$pvalue), ]
resC3_vs_C4Ordered

# Select a subset with padj < 0.05
res_Sub_resC3_vs_C4Ordered <- subset(resC3_vs_C4Ordered, padj < 0.05)
res_Sub_resC3_vs_C4Ordered

# Export the subset with padj < 0.05
write.csv(as.data.frame(res_Sub_resC3_vs_C4Ordered), file = "res_Sub_resC3_vs_C4Ordered[1-0.05].csv")

##############################################################################################

# Data visualization using PCA

# Perform variance stabilizing transformation (VST)
vsd <- vst(ddsTxi.tr, blind = FALSE)
head(assay(vsd), 3)

# Perform logarithm transformation (log2(n + 1))
rld <- rlog(ddsTxi.tr, blind = FALSE)
head(assay(rld), 3)

# Load required libraries
library("dplyr")
library("ggplot2")

# Estimate size factors
dds.tr <- estimateSizeFactors(ddsTxi.tr)

# Prepare data for plotting
df <- bind_rows(
  as_data_frame(log2(counts(dds.tr, normalized = TRUE)[, 1:2] + 1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog")
)

colnames(df)[1:2] <- c("x", "y")
lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels = lvls)

# Plot PCA using ggplot2 and facet by transformation
ggplot(df, aes(x = x, y = y)) +
  geom_hex(bins = 80) +
  coord_fixed() +
  facet_grid(. ~ transformation)

## Calculate sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDists

## Visualize sample distances using a heatmap

library("pheatmap")
library("RColorBrewer")

# Heatmap of sample-to-sample distances using the variance stabilizing transformed values
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Run, vsd$Condition, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)

### PCA plot

plotPCA(vsd, intgroup = c("Sample", "Condition"))

###
pcaData <- plotPCA(vsd, intgroup = c("Sample", "Condition"), returnData = TRUE)
pcaData

##
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition, shape = Sample)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

## PCA plot using Generalized PCA (GPCA) - Uncomment if required

# library("glmpca")
#
# gpca <- glmpca(counts(ddsTxi.tr), L = 2)
# gpca.dat <- gpca$factors
# gpca.dat$Run <- ddsTxi.tr$Run
# gpca.dat$Condition <- ddsTxi.tr$Condition

##################################

# MDS plot

mds <- as.data.frame(colData(vsd)) %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Run, shape = Condition)) +
  geom_point(size = 3) +
  coord_fixed() +
  ggtitle("MDS with VST data")

#######################################################################################

#### Comparison between leaves and cotyledons of Sesuvium sesuviodes

# Load sample file information
setwd("~/RNA_122022")
csvfile1 <- file.path("SampleSSonly.csv")
coldata1 <- read.csv(csvfile1, row.names = 1, stringsAsFactors = FALSE)
coldata1
str(coldata1)

# Convert Condition from character to factor
coldata1$Condition <- as.factor(coldata1$Condition)
str(coldata1)

# Quantification file based on SampleTable information
files1 <- file.path("quant", coldata1$Sample, "quant.sf")
files1

# Assign names to files
names(files1) <- coldata1$Sample
library("readr")

# Import the necessary quantification data for DESeq2 using the tximport function
library("tximport")
library("jsonlite")
txi.tr1 <- tximport(files1, type = "salmon", txOut = TRUE)
txi.tr1
?tximport

# Check count per sample
head(txi.tr$counts)

# Create DESeqDataSet object
ddsTxi.tr1 <- DESeqDataSetFromTximport(txi.tr1, colData = coldata1, design = ~ Condition)

#### Count data (independent of expression analysis) for k-means clustering
# No filter on read count

ddsTxi.tr1 <- DESeq(ddsTxi.tr1)
ddsTxi.tr1

# Retrieve the normalized count matrix
## Normalized count
mat1 <- counts(ddsTxi.tr1, normalized = TRUE)
mat1
write.csv(mat1, file = "norm_countsLeavVScotyNight.csv")

## Unnormalized count
mat11 <- counts(ddsTxi.tr1, normalized = FALSE)
mat11
write.csv(mat11, file = "unnorm_countsLeavVScotyNight.csv")

#################################################################

## Set reference level

ddsTxi.tr1$Condition <- relevel(ddsTxi.tr1$Condition, "CD")

levels(ddsTxi.tr1$Condition)
colData(ddsTxi.tr1)

## Perform differential expression analysis

dds.tr1 <- DESeq(ddsTxi.tr1)
res.tr1 <- results(dds.tr1)
res.tr1
summary(res.tr1)
resultsNames(dds.tr1)

# Extract the general results and order by p-value
resLvsCNigOrdered <- res.tr1[order(res.tr1$pvalue), ]
resLvsCNigOrdered

# Write the general results
write.csv(as.data.frame(resLvsCNigOrdered), file = "resLvsCNigALLOrdered.csv")

## Extract specific results

resCoty_night_vs_Coty_day <- results(dds.tr1, name = "Condition_Coty_night_vs_Coty_day", lfcThreshold = 1, alpha = 0.05)
resCoty_night_vs_Coty_day
summary(resCoty_night_vs_Coty_day)

# Order results by p-value
resCoty_night_vs_Coty_dayOrdered <- resCoty_night_vs_Coty_day[order(resCoty_night_vs_Coty_day$pvalue), ]
resCoty_night_vs_Coty_dayOrdered

# Select a subset with padj < 0.05
res_Sub_resCoty_night_vs_Coty_dayOrdered <- subset(resCoty_night_vs_Coty_dayOrdered, padj < 0.05)
res_Sub_resCoty_night_vs_Coty_dayOrdered

# Export the subset with padj < 0.05
write.csv(as.data.frame(res_Sub_resCoty_night_vs_Coty_dayOrdered), file = "res_Sub_resCoty_night_vs_Coty_dayOrdered[1-0.05].csv")

###############################

# Comparison Condition_Leaf_vs_Coty_day

resLeaf_vs_Coty_day <- results(dds.tr1, name = "Condition_Leaf_vs_Coty_day", lfcThreshold = 1, alpha = 0.05)
resLeaf_vs_Coty_day
summary(resLeaf_vs_Coty_day)

# Order results by p-value
resLeaf_vs_Coty_dayOrdered <- resLeaf_vs_Coty_day[order(resLeaf_vs_Coty_day$pvalue), ]
resLeaf_vs_Coty_dayOrdered

# Select a subset with padj < 0.05
res_Sub_resLeaf_vs_Coty_dayOrdered <- subset(resLeaf_vs_Coty_dayOrdered, padj < 0.05)
res_Sub_resLeaf_vs_Coty_dayOrdered

# Export the subset with padj < 0.05
write.csv(as.data.frame(res_Sub_resLeaf_vs_Coty_dayOrdered), file = "res_Sub_resLeaf_vs_Coty_dayOrdered[1-0.05].csv")

####
# Comparison between Night cotyledons and leaves

## Set night as the reference level

ddsTxi.tr1$Condition <- relevel(ddsTxi.tr1$Condition, "Coty_night")

## Perform differential expression analysis

dds.tr1 <- DESeq(ddsTxi.tr1)
res.tr1 <- results(dds.tr1)
res.tr1
summary(res.tr1)
resultsNames(dds.tr1)

### Comparison night vs leaves

resLeaf_vs_Coty_night <- results(dds.tr1, name = "Condition_Leaf_vs_Coty_night", lfcThreshold = 1, alpha = 0.05)
resLeaf_vs_Coty_night
summary(resLeaf_vs_Coty_night)

# Order results by p-value
resLeaf_vs_Coty_nightOrdered <- resLeaf_vs_Coty_night[order(resLeaf_vs_Coty_night$pvalue), ]
resLeaf_vs_Coty_nightOrdered

# Select a subset with padj < 0.05
res_Sub_resLeaf_vs_Coty_nightOrdered <- subset(resLeaf_vs_Coty_nightOrdered, padj < 0.05)
res_Sub_resLeaf_vs_Coty_nightOrdered

# Export the subset with padj < 0.05
write.csv(as.data.frame(res_Sub_resLeaf_vs_Coty_nightOrdered), file = "res_Sub_resLeaf_vs_Coty_nightOrdered[1-0.05].csv")


#########################################################################
# Data visualization using PCA

# Perform variance stabilizing transformation (VST)
vsd2 <- vst(ddsTxi.tr2, blind = FALSE)
head(assay(vsd2), 3)

###
library("dplyr")
library("ggplot2")

## Samples distance

sampleDists <- dist(t(assay(vsd2)))
sampleDists

## Visualize sample distances using a heatmap

library("pheatmap")
library("RColorBrewer")

# Heatmap of sample-to-sample distances using the variance stabilizing transformed values
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd2$Run, vsd2$Condition, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)

# PCA plot

plotPCA(vsd2, intgroup = c("Sample", "Condition"))

pcaData2 <- plotPCA(vsd2, intgroup = c("Sample", "Condition"), returnData = TRUE)
pcaData2

# Percentage of variance explained by each principal component
percentVar2 <- round(100 * attr(pcaData2, "percentVar"))

# PCA plot with percent variance explained
ggplot(pcaData2, aes(x = PC1, y = PC2, color = Condition, shape = Sample)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar2[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar2[2], "% variance")) +
  coord_fixed() +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
