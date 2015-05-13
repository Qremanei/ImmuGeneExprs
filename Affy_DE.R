# Differential gene expression and gene set enrichment analysis: Microarray
# May 8, 2015
# Qinwen Liu

setwd(base_dir)

########################################
# import microarray ExpressionSet
########################################
library(GEOquery)
gse.eset = getGEO("GSE53166", GSEMatrix=T)
e = gse.eset[[1]]
e.mtx = exprs(e)
str(e.mtx)
sampleInfo = pData(e)
head(sampleInfo)
featureInfo = fData(e)
names(featureInfo)

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
obj.sv = sva(filterEset, mod, mod0)
head(obj.sv$sv)
head(obj.sv$pprob.gam)
obj.sv$n.sv

# perform adjusted differential expression analysis
modSv = cbind(mod, obj.sv$sv)
fit = lmFit(filterEset, modSv)
contrast.matrix <- cbind("C1"=c(0,1,-1, rep(0, 34)),
                         "C2"=c(1,0,-1, rep(0, 34)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
efit <- eBayes(fitContrasts)

tt.1 <- topTable(efit, coef=1, number=dim(filterEset)[1])
idx1 = with(tt.1, {
  logFC > 0.75 | logFC < -1.5 &
  adj.P.Val < 0.01
})
sum(idx1)

tt.2 <- topTable(efit, coef=2, number=dim(filterEset)[1])
idx2 = with(tt.2, {
  logFC > 0.75 | logFC < -1.5 &
  adj.P.Val < 0.01
})
sum(idx2)

idx.comb = unique(c(rownames(tt.1)[idx1], rownames(tt.2)[idx2]))
length(idx.comb)

#####################################################
# Annotating the results with associated gene symbols
#####################################################
library(hugene10sttranscriptcluster.db)
library(annotate)

gene.symbols.1 <- getSYMBOL(rownames(tt.1)[idx1], "hugene10sttranscriptcluster.db")
results1 <- cbind(tt.1[idx1,], gene.symbols.1)
write.table(results1, "LPS_DE_results.txt", sep="\t", quote=FALSE)

gene.symbols.2 <- getSYMBOL(rownames(tt.2)[idx2], "hugene10sttranscriptcluster.db")
results2 <- cbind(tt.2[idx2,], gene.symbols.2)
write.table(results2, "dNS1_DE_results.txt", sep="\t", quote=FALSE)

###################################
# gene set enrichment analysis
###################################
array.db = hugene10sttranscriptcluster.db
columns(array.db)
gene.go = select(array.db, keys=rownames(filterEset), keytypes="PROBEID",
                 columns=c("GO", "SYMBOL"))
head(gene.go)
# length(unique(gene.go$GO))

# enrichment of immune response genes
set.seed(1)
go.idx <- grep("GO:0006955", gene.go$GO)
length(go.idx)
ind = unique(gene.go[go.idx, "PROBEID"])
r1 <- roast(filterEset, ind, design=modSv, contrast=contrast.matrix[,1])
r1
r2 <- roast(filterEset, ind, design=modSv, contrast=contrast.matrix[,2])
r2
# pick a contrast in random as a control
cont.r = c(rep(0, 3), 1, -1, rep(0, 32))
roast(filterEset, ind, design=modSv, contrast=cont.r)

# Testing multiple gene sets
go.group = split(gene.go, gene.go$GO)
head(go.group)
length(go.group)
ind.m = lapply(go.group, function(group) unique(group$PROBEID))
head(ind.m)
length(ind.m)
ind.length = sapply(ind.m, length)
head(ind.length)
ind.go.group = ind.m[ind.length >= 50]
head(ind.go.group)
length(ind.go.group)

set.seed(1)
mr1 <- mroast(filterEset, ind.go.group, design=modSv,  contrast=contrast.matrix[,1])
head(mr1)
mr1 <- mr1[order(-mr1$PropUp),]
head(mr1)

mr2 <- mroast(filterEset, ind.go.group, design=modSv,  contrast=contrast.matrix[,2])
head(mr2)
mr2 <- mr2[order(-mr2$PropUp),]
head(mr2)

# extract the GO terms for the top results, by the mixed test
library(GO.db)
columns(GO.db)
keytypes(GO.db)
GOTERM[[rownames(mr1)[1]]]

mr1tab <- select(GO.db, keys=rownames(mr1)[1:20],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
mr1tab[,1:2]

mr2tab <- select(GO.db, keys=rownames(mr2)[1:20],
                 columns=c("GOID","TERM","DEFINITION"), 
                 keytype="GOID")
mr2tab[,1:2]

mr1tab.up <- select(GO.db, keys=rownames(mr1)[mr1$Direction == "Up"][1:20],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
mr1tab.up[,1:2]

mr1tab.down <- select(GO.db, keys=rownames(mr1)[mr1$Direction == "Down"][1:20],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
mr1tab.down[,1:2]

#####################
sessionInfo()




