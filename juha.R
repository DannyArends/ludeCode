setwd("~/Github/Juha/")

wholeblood <- read.csv("GPL570_WholeBlood.txt",sep='\t',row.names=1)
neutrophil <- read.csv("GPL570_Neutrophil.txt",sep='\t',row.names=1)

#Neutrophil vector lude
ludevectordata <- read.csv("NeutrophilVectorLude.txt",sep="\t")

translation <- read.csv("GPL570ProbeENSGInfo+HGNC.txt",sep='\t',row.names=1)
translation <- translation[which(translation[,9] != "-"),]
ivectorAnn <- read.table("IlluminaProbeTranslationTable.txt",sep='\t',header=TRUE)
ivectorAnn <- ivectorAnn[which(ivectorAnn[,5] != "-"),]

badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")

neutrophil <- neutrophil[,-which(colnames(neutrophil) %in% badneutros)]

IDneutro <- 1:ncol(neutrophil)
IDwhole  <- (ncol(neutrophil)+1):(ncol(wholeblood)+ncol(neutrophil))

combined <- cbind(neutrophil, wholeblood)
combined <- combined[which(rownames(combined) %in% translation[,1]),]

Zscores <- apply(combined, 1, function(x){
  return(
    (mean(x[IDneutro]) - mean(x)) / sd(x)
  )
})

#top1000 <- sort(abs(Zscores))[length(Zscores)-1000]
#Zscores <- Zscores[which(abs(Zscores) > top1000)]

# Add annotation to the cell type vector
sortR <- match(names(Zscores), translation[,1])
myvector <- cbind(translation[sortR,c(1,9)],Zscores) # Add annotation

ludevectordata <- ludevectordata[which(ludevectordata[,1] %in% myvector[,1]),]
sortL <- match(myvector[,1], ludevectordata[,1])

myvector <- cbind(myvector,ludevectordata[sortL,2]) # Add lude's vector

cat("Number significant changed (P < 0.05):", length(which(Zscores < 0.05/length(Zscores))), "\n")

ivector <- read.csv("2013-06-21-EGCUT-Vector-rs12057769-2000128.txt",sep='\t',row.names=NULL)
ivector <- ivector[which(as.character(ivector[,1]) %in% as.character(ivectorAnn[,6])),]

#top1000 <- sort(abs(ivector[,2]))[nrow(ivector)-1000]
#ivector <- ivector[which(abs(ivector[,2]) > top1000),]

# Add annotation to the EGCUT vector
sortA <- match(ivector[,1], ivectorAnn[,6])
ivector <- cbind(ivectorAnn[sortA,5:6], ivector)

# Merge the two vectors
ivector <- ivector[which(ivector[,1] %in% myvector[,2]),]
sortT <- match(ivector[,1],myvector[,2])

resultvector <- cbind(myvector[sortT,c(2,3,4)], ivector[,c(1,4)])

plot(resultvector[,c(2,5)], pch=19, cex=0.4)
cor(resultvector[,c(2,3,5)],method="spearman")

