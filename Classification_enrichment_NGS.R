# Classification and gene set enrichment analysis: NGS
# May 8, 2015
# Qinwen Liu

base_dir = "~/Documents/"
setwd(base_dir)

################################
# import gene expression data
################################

library(GEOquery)
gse.eset = getGEO("GSE53166", GSEMatrix=T)
e = gse.eset[[1]]
e.mtx = exprs(e)
str(e.mtx)
sampleInfo = pData(e)
head(sampleInfo)
featureInfo = fData(e)
head(featureInfo)

#########################################
# EDA
#########################################
library(RColorBrewer)

# quick boxplot: sample a subset of all samples
set.seed(100)
sample_col = sample(1:dim(e.mtx)[2], 30)
sample_row = sample(1:dim(e.mtx)[1], 10000)
cols <- rainbow(12)
boxplot(e.mtx[sample_row, sample_col], col=cols)

# MA plot of a subset of the probe intensity measurement
library(rafalib)
mypar(2,5)
sapply(sample_col[2:11], function(n) {
  M = e.mtx[sample_row, n] - e.mtx[sample_row, sample_col[1]]
  A = (e.mtx[sample_row, n] + e.mtx[sample_row, sample_col[1]])/2
  plot(A, M, xlab=paste("A:sample_1", n, sep="_"), ylab=paste("M:sample_1", n, sep="_"), pch=15, cex=0.5)
  abline(h=0)
})
dev.off()

###############################################################
# filter data: remove genes with little cross-sample variation
###############################################################
library(genefilter)
e.var = rowVars(e.mtx)
hist(e.var, breaks=100)

e.filtered <- nsFilter(e, require.entrez=FALSE, remove.dupEntrez=FALSE)
# What got removed and why
e.filtered$filter.log
filterEset <- exprs(e.filtered$eset)
dim(filterEset)

library(MLSeq)
library(DESeq2)

filepath = system.file("extdata/cervical.txt", package = "MLSeq")
cervical = read.table(filepath, header = TRUE)
str(cervical)

set.seed(9)

class = data.frame(condition = factor(rep(c(0, 1), c(29, 29))))

# select test set and its sample class
nTest = ceiling(ncol(cervical) * 0.2)
ind = sample(ncol(cervical), nTest, FALSE)
cervical.test = cervical[, ind]
cervical.test = as.matrix(cervical.test + 1)
classts = data.frame(condition = class[ind, ])

# select training set and its sample class
cervical.train = cervical[, -ind]
cervical.train = as.matrix(cervical.train + 1)
classtr = data.frame(condition = class[-ind, ])

cervical.trainS4 = DESeqDataSetFromMatrix(countData = cervical.train, 
                                          colData = classtr, formula(~condition))
class(cervical.trainS4)
cervical.trainS4 = DESeq(cervical.trainS4, fitType = "local")

cervical.testS4 = DESeqDataSetFromMatrix(countData = cervical.test, colData = classts,
                                         formula(~condition))
cervical.testS4 = DESeq(cervical.testS4, fitType = "local")

# Classify using Support Vector Machines
svm = classify(data = cervical.trainS4, method = "svm", normalize = "deseq",
               deseqTransform = "vst", cv = 5, rpt = 3, ref = "1")
svm
getSlots("MLSeq")
trained(svm)

# predict the class labels
pred.svm = predictClassify(svm, cervical.testS4)
table(pred.svm, relevel(cervical.testS4$condition, 2))

############################
# gene set analysis
############################

library(goseq)
library(edgeR)
path <- system.file(package="goseq", "extdata", "Li_sum.txt")

table.summary <- read.table(path, sep='\t', header=TRUE, stringsAsFactors=FALSE)
counts <- table.summary[,-1]
rownames(counts) <- table.summary[,1]
grp <- factor(rep(c("Control","Treated"), times=c(4,3)))
summarized <- DGEList(counts, lib.size=colSums(counts), group=grp)

# Use a ‘common’ dispersion estimate, and compare the two groups using an exact test
disp <- estimateCommonDisp(summarized)
tested <- exactTest(disp)
topTags(tested)

padj <- with(tested$table, {
  keep <- logFC != 0
  value <- p.adjust(PValue[keep], method="BH")
  setNames(value, rownames(tested)[keep])
})
genes <- padj < 0.05
table(genes)

pwf <- nullp(genes,"hg19","ensGene")
head(pwf)

# association of genes to GO pathway
GO.wall <- goseq(pwf, "hg19", "ensGene")
head(GO.wall)

GO.nobias <- goseq(pwf,"hg19","ensGene",method="Hypergeometric")

#Compare the over-represented P-values for each set, under the different methods
idx <- match(GO.nobias$category, GO.wall$category)
plot(log10(GO.nobias[, "over_represented_pvalue"]) ~
       log10(GO.wall[idx, "over_represented_pvalue"]),
     xlab="Wallenius", ylab="Hypergeometric",
     xlim=c(-5, 0), ylim=c(-5, 0))
abline(0, 1, col="red", lwd=2)

####################
# limma package
####################

library(GEOquery)
g <- getGEO("GSE34313")
e <- g[[1]]

e$condition <- e$characteristics_ch1.2
levels(e$condition) <- c("dex24","dex4","control")
table(e$condition)
boxplot(exprs(e), range=0)

names(fData(e))
lvls <- c("control", "dex4")
es <- e[,e$condition %in% lvls]
es$condition <- factor(es$condition, levels=lvls)

library(limma)
design <- model.matrix(~ es$condition)
fit <- lmFit(es, design=design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, genelist=fData(es)$GENE_SYMBOL)
tt
topTable(fit)[,c(6,7,18,22)]

set.seed(1)
idx <- grep("GO:0045454", fData(es)$GO_ID)
length(idx)
r1 <- roast(es, idx, design)
?roast
r1

# Testing multiple gene sets
library(org.Hs.eg.db)
org.Hs.egGO2EG
go2eg <- as.list(org.Hs.egGO2EG)
head(go2eg)

govector <- unlist(go2eg)
golengths <- sapply(go2eg, length)
head(fData(es)$GENE)
# map entrez gene ID to row numbers of expression set
idxvector <- match(govector, fData(es)$GENE)
table(is.na(idxvector))
# row index in expression set
idx <- split(idxvector, rep(names(go2eg), golengths))
go2eg[[1]]
fData(es)$GENE[idx[[1]]]

idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths >= 50]
length(idxsub)

set.seed(1)
r2 <- mroast(es, idxsub, design)
head(r2)
r2 <- r2[order(-r2$PropUp),]
head(r2)

# extract the GO terms for the top results, by the mixed test
library(GO.db)
columns(GO.db)
keytypes(GO.db)
GOTERM[[rownames(r2)[1]]]

select(GO.db, keys=rownames(r2)[1], columns="TERM", keytype="GOID")

r2tab <- select(GO.db, keys=rownames(r2)[1:10],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]

r2 <- r2[order(r2$PValue),]
r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction == "Up"][1:10],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]

r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction == "Down"][1:5],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]

library(MASS)
Sigma = matrix(.7, ncol=10, nrow=10)
diag(Sigma) = 1
mtx = mvrnorm(n=10000,mu=rep(0,10),Sigma=Sigma)
var(rowMeans(mtx))






