setwd("~/Github/Juha/")

#metaRes <- read.csv("MetaAnalysisZScoreMatrix-Ensembl.txt",sep='\t',row.names=1)
load("metaRes.Rdata")

wholeblood <- read.csv("GPL570_WholeBlood.txt", sep='\t', row.names=1)
Neutr      <- read.csv("GPL570_Neutrophil.txt", sep='\t', row.names=1)
Bcell      <- read.csv("GPL570_Bcell.txt", sep='\t', row.names=1)
Tcell      <- read.csv("GPL570_Tcell.txt", sep='\t', row.names=1)
NKcell     <- read.csv("GPL570_NKcell.txt", sep='\t', row.names=1)
RBC        <- read.csv("GPL570_RBC.txt", sep='\t', row.names=1)

translation <- read.csv("GPL570ProbeENSGInfo+HGNC.txt",sep='\t',row.names=1)
translation <- translation[which(translation[,9] != "-"),]
ivectorAnn <- read.table("IlluminaProbeTranslationTable.txt",sep='\t',header=TRUE)
ivectorAnn <- ivectorAnn[which(ivectorAnn[,5] != "-"),]

LudeVec <- read.csv("NeutrophilVectorLude.txt",sep="\t")
IntrVec <- read.csv("2013-06-21-EGCUT-Vector-rs12057769-2000128.txt",sep='\t',row.names=NULL)

badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")

Neutr  <- Neutr[,-which(colnames(Neutr) %in% badneutros)]  # Remove remaining bad neutrophil samples
Bcell  <- Bcell[,-which(apply(cor(Bcell), 2, median) < 0.9)] # Use correlation for Bcell
Tcell  <- Tcell[,-which(apply(cor(Tcell), 2, median) < 0.9)] # Use correlation for Tcell

IntrVec <- IntrVec[which(as.character(IntrVec[,1]) %in% as.character(ivectorAnn[,6])),]
row <- 0
zscores <- t(apply(wholeblood, 1, function(x){ # Using basis T statistics actually
  row <<- row + 1
  meanWB <- mean(x)
  sdWB <- sd(x)
  if(row %% 100 == 0)cat("Done:", 3*row,"t tests\n")
  return( cbind(t.test(as.numeric(Neutr[row,]), x)$statistic, 
                t.test(as.numeric(Bcell[row,]), x)$statistic, 
                t.test(as.numeric(Tcell[row,]), x)$statistic,
                t.test(as.numeric(NKcell[row,]),x)$statistic,
                t.test(as.numeric(RBC[row,]),x)$statistic) )
}))
colnames(zscores) <- c("Neutrophil", "Bcell", "Tcell", "NKcell", "RBC")

# Add annotation to the cell type vector
zscores <- zscores[which(rownames(zscores) %in% translation[,1]),]
sortZ <- match(rownames(zscores), translation[,1])
scores <- cbind(translation[sortZ, c(1,9)], zscores) # Add annotation

# Add annotation to the EGCUT vector
sortA <- match(IntrVec[,1], ivectorAnn[,6])
IntrVec <- cbind(ivectorAnn[sortA, 5:6], IntrVec)

# Merge the two vectors
IntrVec <- IntrVec[which(IntrVec[,1] %in% scores[,2]),]
sortT <- match(IntrVec[,1], scores[,2])

ResVector <- cbind(scores[sortT,], IntrVec[,c(1,4)])

# Bind the ResVector columns Neutrophil, Bcell and Tcell to the MetaAnalysis data
metaRes <- metaRes[which(rownames(metaRes) %in% ResVector[,1]),]     # match MetaRes
ResVector <- ResVector[which(ResVector[,1] %in% rownames(metaRes)),] # match ResVector
sortRes <- match(rownames(metaRes), ResVector[,1]) # Align

metaRes <- cbind(ResVector[sortRes, 3:7], metaRes)

corrs <- NULL  # Correlation analysis between 1:3 Neutro, B, Tcell, and all other columns
cnt <- 1
snpnames <- colnames(metaRes[, 6:length(metaRes)])
r <- apply(metaRes[, 6:length(metaRes)], 2, function(qtl){
  corrs <<- rbind(corrs, cor(qtl, metaRes[,1:5],use="pair", method="spearman"))
  rownames(corrs) <<- snpnames[1:cnt]
  cnt <<- cnt + 1
  if(cnt %% 100 == 0)cat("Cor of", cnt,"(QTL:cType):SNP pairs\n")
})

write.csv(corrs, file="COR_MetaAnalysisZScore_CellTypeTScore.txt", quote = FALSE)

summary(corrs[,1]) # Neutrophil -0.6 to 0.6
summary(corrs[,2]) # Bcell      -0.4 to 0.4
summary(corrs[,3]) # Tcell      -0.6 to 0.6
summary(corrs[,4]) # NKcell     -0.5 to 0.5
summary(corrs[,5]) # RBC        -0.2 to 0.2

top100_Neutro <- rownames(corrs)[sort(abs(corrs[,1]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_Bcell <- rownames(corrs)[sort(abs(corrs[,2]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_Tcell <- rownames(corrs)[sort(abs(corrs[,3]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_NKcell <- rownames(corrs)[sort(abs(corrs[,4]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]
top100_RBC <- rownames(corrs)[sort(abs(corrs[,5]),index.ret=TRUE, decreasing=TRUE)$ix[1:100]]

heatmap(corrs[c(top100_Neutro,top100_Bcell,top100_Tcell,top100_NKcell,top100_RBC), ], col=c("red","white","green"))

top <- which(apply(corrs,1,function(x){max(abs(x)) > 0.4}))
aa <- heatmap(corrs[top,],keep.dendro=TRUE)
heatmap.2(t(corrs[top[aa$rowInd],]), trace="none", col=c("red","pink","white","lightgreen","green"), main="Correlation CellTypes vs QTL Z-scores", key=FALSE, margins=c(5,10), scale="none", dendrogram="row", labCol="",xlab="Zscore Snp:Probe")

#plot(ResVector[,c(2,5)], pch=19, cex=0.4)
cor(ResVector[,c(3, 4, 5, 7)], use="pair",method="spearman")

write.csv(ResVector[,-6], file="Tvector.txt", quote = FALSE)

