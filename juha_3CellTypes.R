setwd("~/Github/Juha/")

metaRes <- read.csv("MetaAnalysisZScoreMatrix-Ensembl.txt",sep='\t',row.names=1)

wholeblood <- read.csv("GPL570_WholeBlood.txt",sep='\t',row.names=1)
Neutr <- read.csv("GPL570_Neutrophil.txt",sep='\t',row.names=1)
Bcell <- read.csv("GPL570_Bcell.txt",sep='\t',row.names=1)
Tcell <- read.csv("GPL570_Tcell.txt",sep='\t',row.names=1)

translation <- read.csv("GPL570ProbeENSGInfo+HGNC.txt",sep='\t',row.names=1)
translation <- translation[which(translation[,9] != "-"),]
ivectorAnn <- read.table("IlluminaProbeTranslationTable.txt",sep='\t',header=TRUE)
ivectorAnn <- ivectorAnn[which(ivectorAnn[,5] != "-"),]

LudeVec <- read.csv("NeutrophilVectorLude.txt",sep="\t")
IntrVec <- read.csv("2013-06-21-EGCUT-Vector-rs12057769-2000128.txt",sep='\t',row.names=NULL)


badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")

neutrophil <- neutrophil[,-which(colnames(neutrophil) %in% badneutros)]

Bcell <- Bcell[,-which(apply(cor(Bcell),2,median) < 0.9)]
Tcell <- Tcell[,-which(apply(cor(Tcell),2,median) < 0.9)]

IntrVec <- IntrVec[which(as.character(IntrVec[,1]) %in% as.character(ivectorAnn[,6])),]
row <- 0
zscores <- t(apply(wholeblood, 1, function(x){
  row <<- row+1
  meanWB <- mean(x)
  sdWB <- sd(x)
  return( cbind(t.test(as.numeric(Neutr[row,]), x)$statistic, 
                t.test(as.numeric(Bcell[row,]), x)$statistic, 
                t.test(as.numeric(Tcell[row,]), x)$statistic))
}))
colnames(zscores) <- c("Neutrophil", "Bcell", "Tcell")

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

#plot(ResVector[,c(2,5)], pch=19, cex=0.4)
cor(ResVector[,c(3, 4, 5, 7)], use="pair",method="spearman")

write.csv(ResVector[,-6], file="Tvector.txt", quote = FALSE)

