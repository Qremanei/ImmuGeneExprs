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


