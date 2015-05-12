# Affymetrix microarray data analysis
# March 31, 2015

setwd(base_dir)

#####################################
# install and load required packages
#####################################

# package for getting GEO datasets
library(GEOquery)
# package for analyzing Affymetrix microarray
library(simpleaffy)

library(RColorBrewer)
library(affyPLM)

# install rgl package for PCA plotting
source("http://bioconductor.org/biocLite.R")
biocLite("rgl")
library(rgl)

library(limma)

# install annotation packages
# the particular Affymetrix chip: hgu133plus2.db
biocLite("hgu133plus2.db")
library(hgu133plus2.db)
library(annotate)

#####################################
# get the data and 
#####################################
#library(GEOquery)
getGEOSuppFiles("GSE20986")
untar("GSE20986/GSE20986_RAW.tar", exdir="GSE20986/data")
cels <- list.files("GSE20986/data/", pattern= "[gz]")
cels
# CEL files (the Affymetrix native file format)
sapply(paste("GSE20986/data", cels, sep="/"), gunzip)
cels

#####################################
# describe the experiment
#####################################

# create a txt file containing Name, FileName, Target

#####################################
# load and normalize the data
#####################################
#library(simpleaffy)
setwd(dir)
celfiles <- read.affy(covdesc="phenodata.txt", path="GSE20986/data")
celfiles

# GC-RMA algorithm is a good method for Affymetrix chip by experience
celfiles.gcrma <- gcrma(celfiles)

celfiles.gcrma

#####################################
# quality control checks:
#####################################
# library(RColorBrewer)
# plot a boxplot of unnormalized intensity values
cols <- rainbow(12)
boxplot(celfiles, col=cols)
hist(celfiles, col=cols)

# plot a boxplot of normalized intensity values
# library(affyPLM)
boxplot(celfiles.gcrma, col=cols)
hist(celfiles.gcrma, col=cols)

eset <- exprs(celfiles.gcrma)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters)

#####################################
# filter data
#####################################
library(genefilter)
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE,
                              remove.dupEntrez=FALSE)
# What got removed and why
celfiles.filtered$filter.log
filterEset <- exprs(celfiles.filtered$eset)
dim(filterEset)

#####################################
# Exploratory analysis by PCA
#####################################
# library(rgl)
pca <- prcomp(t(filterEset), scale=TRUE)
summary(pca)

myColors <- c("Blue", "yellow", "yellow", "Blue", "yellow", "Blue",
              rep("Red", 3), rep("Green", 3))
plot3d(pca$x[, 1:3], col=myColors, xlab="PC1", ylab = "PC2", zlab =
         "PC3", type = "s")

# graphic device has to be open for the code to work
rgl.postscript("PCA.pdf", fmt="pdf", drawText=TRUE)

samples <- celfiles.gcrma$Target
samples <- as.factor(samples)

# set up the experimental design
design <- model.matrix(~0 + samples)
colnames(design) <- c("choroid", "huvec", "iris", "retina")
design

# library(limma)
# fit the linear model to the filtered expression set
fit <- lmFit(filterEset, design)
contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid,
                                 huvec_retina = huvec - retina, huvec_iris <- huvec - iris,
                                 levels=design)

# Now the contrast matrix is combined with the per-probeset linear model fit
huvec_fits <- contrasts.fit(fit, contrast.matrix)

huvec_ebFit <- eBayes(huvec_fits)
# coef=1 is huvec_choroid, coef=2 is huvec_retina, coef=3 is huvec_iris
topTable(huvec_ebFit, number=10, coef=1)
topTable(huvec_ebFit, number=20, coef=3)
topTable(huvec_ebFit, coef=1, number=10000, lfc=5)
nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=5))

probeset.list <- topTable(huvec_ebFit, coef=1, p.value=0.05, lfc=5)
selSamples <- subset(celfiles.gcrma$FileName, (celfiles.gcrma$Target
                                               == "choroid") | (celfiles.gcrma$Target == "huvec"))
probeset.list1 <- topTable(huvec_ebFit, coef=1,
                           number=nrow(filterEset), lfc=0)
probeset.list2 <- probeset.list1[(probeset.list1$adj.P.Val <= 0.05) &
                                   (abs(probeset.list1$logFC) >= 2), ]

selData <- filterEset[rownames(filterEset) %in%
                        rownames(probeset.list2), colnames(filterEset) %in% selSamples]
pdf(file="Heatmap.pdf", width=8, height=6)
heatmap(selData, labRow=c(""), col=topo.colors(16), cexCol=0.6)
graphics.off()

#####################################################
# Annotating the results with associated gene symbols
#####################################################
# library(hgu133plus2.db)
# library(annotate)
gene.symbols <- getSYMBOL(rownames(probeset.list), "hgu133plus2")
results <- cbind(probeset.list, gene.symbols)
write.table(results, "results.txt", sep="\t", quote=FALSE)

save.image()
# Print version information about R, the OS and attached or loaded packages
sessionInfo()






















