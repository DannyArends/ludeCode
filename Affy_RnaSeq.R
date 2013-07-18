setwd("~/Github/Juha/")

RNASeq <- read.csv("expression_table.genes.exonic_v69.0.3.rawCounts.txt", sep='\t', row.names=1)
sampleNames <- read.csv("SampleDetails.txt", sep='\t', row.names=1)
colnames(RNASeq) <- sampleNames[,1]

ordering <- unlist(lapply(strsplit(colnames(RNASeq),"_"),"[",6))

# Neutrophils
Neutr      <- read.csv("GPL570_Neutrophil.txt", sep='\t', row.names=1)

badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")

Neutr  <- Neutr[,-which(colnames(Neutr) %in% badneutros)]  # Remove remaining bad neutrophil samples

# Translation
translation <- read.csv("GPL570ProbeENSGInfo+HGNC.txt",sep='\t',row.names=1)
translation <- translation[which(translation[,9] != "-"),]

# Translation
inTrans <- which(rownames(Neutr) %in% translation[,1])
Neutr <- Neutr[inTrans,]
sortTrans <- match(rownames(Neutr), translation[,1]) # Align
Neutr <- cbind(translation[sortTrans,9], Neutr)

inSeq   <- which(as.character(Neutr[,1]) %in% rownames(RNASeq))
Neutr   <- Neutr[inSeq,]
sortSeq <- match(as.character(Neutr[,1]), rownames(RNASeq)) # Align

Neutr <- cbind(RNASeq[sortSeq,], Neutr)
write.table(Neutr,"Neutrophil_RNAseq.txt",sep='\t',quote=FALSE)
NeutrRna <- log2(Neutr[,7])
keep <- which(is.finite(NeutrRna)) #Kepp onlu the finite ones
Neutr <- Neutr[keep, ]
NeutrRna <- NeutrRna[keep]
NeutrMed <- apply(Neutr[,9:ncol(Neutr)], 1, median)
NeutrMea <- apply(Neutr[,9:ncol(Neutr)], 1, mean)

lowexpressed <- which(NeutrMea < 4 | NeutrRna < 3)
NeutrRna <- NeutrRna[-lowexpressed]
NeutrMea <- NeutrMea[-lowexpressed]
NeutrMed <- NeutrMed[-lowexpressed]
Neutr <- Neutr[-lowexpressed,]

cmea <- cor(NeutrMea, NeutrRna, method="spearman")
cmed <- cor(NeutrMed, NeutrRna, method="spearman")
colmea <- as.numeric(abs(NeutrRna/max(NeutrRna) - NeutrMea/max(NeutrMea)) > 0.3)+1
colmed <- as.numeric(abs(NeutrRna/max(NeutrRna) - NeutrMed/max(NeutrMed)) > 0.3)+1
plot(NeutrMea, NeutrRna, xlab = "Affy", ylab = "RNAseq",main=paste0("Mean Cor: ",cmea), col = colmea, cex=0.7)
plot(NeutrMed, NeutrRna, xlab = "Affy", ylab = "RNAseq",main=paste0("Median Cor: ",cmed), col = colmed, cex=0.7)

# OLD: c(which(NeutrMea > 10 & NeutrRna < 10), which(NeutrMea > 8 & NeutrRna < 5))
belowBanana <- which(colmea==2)
below <- Neutr[belowBanana,]

plot(c(0,nrow(below)), c(0,20),t='n')
boxplot(t(below[,9:ncol(below)]),add=TRUE)
points(log2(below[,7]), col='red')
dev.off()

cat(names(belowBanana), sep='\n')

