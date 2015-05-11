# Gene expression modeling and gene set enrichment analysis
# May 8, 2015
# Qinwen Liu

setwd(base_dir)

#################################################
# import microarray or NGS gene expression data
#################################################

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

#####################################
# Exploratory analysis by PCA
#####################################
library(rgl)

pca <- prcomp(t(filterEset), scale=TRUE)
summary(pca)

myColors <- as.numeric(sampleInfo$characteristics_ch1)
plot3d(pca$x[, 1:3], col=myColors, xlab="PC1", ylab = "PC2", zlab =
         "PC3", type = "s")

# graphic device has to be open for the code to work
# rgl.postscript("PCA.pdf", fmt="pdf", drawText=TRUE)

############################################################
# surrogate variable analysis and regression modeling
############################################################
library(limma)
library(sva)

head(sampleInfo)
# create the full model matrix containing variables of interest and adjustment variables
mod = model.matrix(~  0 + characteristics_ch1 + characteristics_ch1.5 + characteristics_ch1.9 +
                     characteristics_ch1.10, data=sampleInfo)
# create the null model matrix containing adjustment variable
mod0 = model.matrix(~ 0 + characteristics_ch1.5 + characteristics_ch1.9 + characteristics_ch1.10 - 1,
                    data=sampleInfo)

# estimate the surrogate variables
obj.sv = sva(e.mtx, mod, mod0)
head(obj.sv$sv)
head(obj.sv$pprob.gam)
obj.sv$n.sv

# perform adjusted differential expression analysis
modSv = cbind(mod, obj.sv$sv)
fit = lmFit(e.mtx, modSv)
contrast.matrix <- cbind("C1"=c(0,1,-1, rep(0, 34)),
                         "C2"=c(1,0,-1, rep(0, 34)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
efit <- eBayes(fitContrasts)

tt <- topTable(efit, adjust="BH")
tt


###################################
# gene set enrichment analysis
###################################

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


sessionInfo()




