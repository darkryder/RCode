require(xtable)
source("Functions.R")

###################################################
### LoadData
###################################################
resDir <- ".Eh_"

### Read raw counts data
inputFile <- "Eh.txt"
rawdata <- read.table(file.path(inputFile),header=TRUE,sep="\t")
inputFilePrint <- gsub("_","-",inputFile)

## Data normalization??
geneCount <- rawdata[,-c(1,8)]
rc <- rowSums(geneCount)
geneCount <- geneCount[rc > 0,]
rownames(geneCount) <- rawdata[rc > 0,1]
N <- colSums(geneCount)

##Trs Length
trsLength <- NA
if(is.element("length", colnames(rawdata)))
    trsLength <- rawdata[rc > 0,"length"]

##Group
group <- factor(c("A","B","A","A","B","B"))
condA="A"
condB="B"

###################################################
### Normalization
###################################################??
norm <- getNormData(geneCount, group, trsLength, addRaw=TRUE)
save(norm, file=paste(resDir, "_norm.RData", sep=""))
normBoxplot(norm)

###################################################
### VarIntra
###################################################
varIntra(norm, group)

###################################################
### Differential Analysis
###################################################
result <- diffAnDESseq(geneCount, group, trsLength, round=FALSE, condA, condB)
DEG <- lapply(result, function(r){r[which(r$padj<=0.05),"id"]})
nbDEG <- unlist(lapply(DEG, length))

## calcul pour matrice consensus
mat <- matrix(0, nrow=nrow(result[[1]]), ncol=length(result))
for (i in 1:length(result)) mat[,i] <- result[[i]][,"padj"]
colnames(mat) <- names(result)
rownames(mat) <- result[[1]][,"id"]
save(mat, file=paste(resDir, "_matrix.RData", sep=""))

##Table of nbDEG
print(xtable(t(matrix(nbDEG, dimnames=list(names(DEG),"nbDEG"))), caption="Total number of differentially expressed genes for each normalization method (the DESeq package is used for the statistical test of differential expression)", align=c("|c|",rep("c|",length(nbDEG)))), "latex", size="small")

###################################################
### Comparison between list of DEG
###################################################
tabcomp <- matrix(NA, ncol=length(DEG), nrow=length(DEG))
rownames(tabcomp) <- colnames(tabcomp)<-names(DEG)
for (i in names(DEG)) for (j in names(DEG)) tabcomp[i,j]<-length(intersect(DEG[[i]],DEG[[j]]))
print(xtable(tabcomp, caption="Differentially expressed gene list comparison", align=c("|c|",rep("c|",ncol(tabcomp)))), "latex", size="small")

###################################################
### Hierarchical Clustering
###################################################
hclustDE(DEG)

###################################################
### MAPlot
###################################################
MAplotDE(geneCount, group, condA, condB, DEG , showComm=TRUE, showSpe=TRUE)

