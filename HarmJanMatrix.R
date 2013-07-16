setwd("~/Github/Juha/")

ExpCorMatrixHJ <- read.csv("EGCUTEndophenotypesValidSamples-NoSex.txt-asLoadedByNormalizerEndophenotypeVsExpressionCorrelationMatrix-Ensembl.txt", sep='\t', row.names=1)
TScores <- read.csv("TScores5CellTypes.txt", row.names=1)

ExpCorMatrixHJ <- ExpCorMatrixHJ[which(rownames(ExpCorMatrixHJ) %in% rownames(Tscores)),]
TScores <- TScores[which(rownames(TScores) %in% rownames(ExpCorMatrixHJ)),]

sortHJ <- match(rownames(TScores), rownames(ExpCorMatrixHJ))
bound <- cbind(ExpCorMatrixHJ[sortHJ,], TScores)

corrs <- cor(bound, method="spearman")
rownames(corrs)[20:24] <- paste0(rownames(corrs)[20:24],"_Danny")
write.table(corrs,file="Tscores_vs_EGCUTEndophenotypesValidSamples.txt",sep='\t',quote=FALSE)

namez <- colnames(bound)
for(x in 20:24){
  for(y in 10:14){
    xlab=paste0("Tscore ",namez[x])
    ylab=paste0("Endophenotype ",namez[y])
    main=paste0(namez[x]," ", namez[y], ", Cor: ", round(corrs[x,y],d=2))
    png(file=paste0(namez[x],"_",namez[y],".png"),width=1024,height=800)
    plot(bound[,x], bound[,y], xlab=xlab, ylab=ylab, main=main, pch=19, cex=0.6)
    dev.off()
  }
}

