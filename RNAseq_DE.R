# pipeline for analyzing RNAseq datasets
# Differential gene expression analysis: RNAseq
# May 4, 2015
# Qinwen Liu

###############################
# import alignment bam files
###############################

library(Rsamtools)

# set up directories for bam files and metadata table
dir_bamfiles <- "~/Documents/trainings/bioinfo2014/RNAseq/RNAseq/cri2014_rnaseq/bamfiles"
list.files(dir_bamfiles)

dir_sampleTable <- "~/Documents/trainings/bioinfo2014/RNAseq/RNAseq/sampleinfo.txt"
sampleTable <- read.table(dir_sampleTable, header=T, sep="\t")
sampleTable
filenames = file.path(dir_bamfiles, paste0(sampleTable$Sample, ".tophat.bam"))
filenames

bamfiles <- BamFileList(filenames, yieldSize=2000000)
# check the chromosome names in the alignment files
seqinfo(bamfiles)

# read in gene models, and convert GTF file to a TranscriptDb object (TxDb)
library(GenomicFeatures)
gtffile <- "~/Documents/trainings/bioinfo2014/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.gtf"
txdb <- makeTranscriptDbFromGFF(gtffile, format="gtf")
txdb
genes <- exonsBy(txdb, by="gene")
genes

# bioconductor TxDB package
# library(BiocInstaller)
# biocLite("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
# library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
# TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

##########################
# RNAseq read counting
##########################

# count the number of reads unambiguously assigned to genomic features for each sample
library(GenomicAlignments)

se <- summarizeOverlaps(features=genes, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
se
# access the SummarizedExperiment class
dim(se)
dimnames(se)
assay(se)[1:5, 1:6]

rowData(se)
metadata(rowData(se))

colData(se)
colData(se) <- DataFrame(sampleTable, row.names = sampleTable$Sample)
colData(se)
rownames(colData(se))

exptData(se)

#########################################
# Read count EDA
#########################################

round( colSums(assay(se)) / 1e6, 1 )

# biocLite("DESeq2")
library(DESeq2)

dds <- DESeqDataSet(se, design = ~ Group)
dds

# regularized-log transformation
rld <- rlog(dds)
rld
head(assay(rld))
metadata(rowData(rld))

# basic exploratory data analysis of the counts
par( mfrow = c( 1, 2 ) )
?estimateSizeFactors
ddsSF <- estimateSizeFactors(dds)
sizeFactors(ddsSF)
plot( log2( 1 + counts(ddsSF, normalized=TRUE)[ , 1:2] ),
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )

# assess overall similarity between samples: Euclidean distance
sampleDists <- dist( t( assay(rld) ) )
sampleDists

# heatmap
library("gplots")
library("RColorBrewer")

# calculating sample distances is to use the Poisson Distance
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$Group
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

# calculating sample distances is to use the Poisson Distance
install.packages("PoiClaClu")
library(PoiClaClu)

poisd <- PoissonDistance(t(counts(dds)))
names(poisd)
poisd$dd
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- dds$Group
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(poisd$dd)
heatmap.2( samplePoisDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

# PCA: to see NGS batch effect
names(colData(rld))
plotPCA(rld, intgroup = "Group")
data <- plotPCA(rld, intgroup = "Group", returnData=TRUE)

# multidimensional scaling (MDS): use distance matrix
library(ggplot2)
mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, as.data.frame(colData(rld)))
qplot(X1,X2,color=Group,data=mds)

mds <- data.frame(cmdscale(samplePoisDistMatrix))
mds <- cbind(mds, as.data.frame(colData(dds)))
qplot(X1,X2,color=Group,data=mds)

####################################
# Differential expression analysis
####################################

# set the reference level for fold change comparison
dds$Group = relevel(dds$Group, "batch")
dds$Group

## Differential expression analysis based on the Negative Binomial distribution
DEs <- DESeq(dds)
# results in a DESeqDataSet object
DEs
# extract out results tables: a DataFrame object
res <- results(DEs)
mcols(res, use.names=TRUE)
summary(res)
# correct for Multiple testing
resSig <- subset(res, padj < 0.1)
head(resSig[ order( resSig$log2FoldChange ), ])
head(resSig[ order( resSig$log2FoldChange, decreasing=T ), ])

# Hypothetical: add contrast if there are other factors that may affect read counts
# results(DEs, contrast=c("SequencingCenter", "CUT1", "CUT2"))

####################################
# Diagnostic plots
####################################

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(DEs, gene=topGene, intgroup=c("Group"))

# nicer plot
data <- plotCounts(DEs, gene=topGene, intgroup=c("Group"), returnData=TRUE)
ggplot(data, aes(x=Group, y=count, fill=Group)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")

# MA plot: DESeq2 incorporates a prior on log2 fold changes; red: padj < 0.1
plotMA(res, ylim=c(-10,10))
# label the most significant change
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# Whether a gene is called significant depends not only on its LFC (log2 fold change),
# but also on its within-group variability, which DESeq2 quantifies as the dispersion
plotDispEsts(DEs)

# histogram
hist(res$pvalue, breaks=20, col="grey50", border="white")
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")

# Gene clustering
library("genefilter")
topVarGenes <- head(order(-rowVars(assay(rld))),35)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Group ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- rld$Group
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")

# for RNA-seq data, weakly expressed genes have no chance of showing DE because of high Poisson noise
# create bins using the quantile function
qs <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:7/7))
# cut the genes into the bins
bins <- cut(res$baseMean, qs)
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of $p$ values less than .01 for each bin
ratios <- tapply(res$pvalue, bins, function(p) mean(p < .01, na.rm=TRUE))
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")

## Independent filtering:
# By removing the weakly-expressed genes from the input to the FDR procedure
# more genes to be significant among those which we keep so improved the power of the test

attr(res,"filterThreshold")
attr(res, "filterNumRej")
plot(attr(res,"filterNumRej"),type="b",
     xlab="quantiles of 'baseMean'",
     ylab="number of rejections")

####################################
# gene annotation
####################################
biocLite("org.Sc.sgd.db")
library("AnnotationDbi")
library(org.Sc.sgd.db)
org.Sc.sgd.db

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

columns(org.Sc.sgd.db)
res$hgnc_symbol <- convertIDs(row.names(res), "ENSEMBL", "GENENAME", org.Sc.sgd.db)
res$entrezgene <- convertIDs(row.names(res), "ENSEMBL", "ENTREZID", org.Sc.sgd.db)

resOrdered <- res[order(res$pvalue),]
head(resOrdered)

write.csv(as.data.frame(resOrdered), file="results.csv")



















