#Match up
setwd("~/Github/LudeNew")
library(biomaRt)
#mart   <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Use annotation to select good probes #
annot <- read.table("2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",sep='\t',header=TRUE)

removeMissing <- function(matrix, col="HUGO"){
  matrix <- matrix[which(matrix[,col] != "-"),]
  matrix <- matrix[which(matrix[,col] != '-1'), ]
  matrix <- matrix[which(matrix[,col] != ''), ] #Remove Non Annotated probes  
  invisible(matrix)
}

annot <- removeMissing(annot)
#Load in the HarmJan data (atm we only have 1 cell type there
CTVHJ <- NULL ; colHeaders <- NULL
for(x in dir("BloedCelDataHarmJan/")[4]){
  fn <- paste0("BloedCelDataHarmJan/",x)
  rowH <- read.csv(fn, sep='\t')[,1:2]
  colHeaders <- c(colHeaders, x)
  CTVHJ <- cbind(CTVHJ, read.csv(fn, sep='\t')[,4])
}
rownames(CTVHJ) <- rowH[,1]
CTVHJ <- cbind(rowH[,2], CTVHJ)
colnames(CTVHJ) <- c("HT12v3Probe",colHeaders)

# Throw away annotation we dont have a HJ probe for
matchAnnot <-which(as.character(annot[,6]) %in% as.character(CTVHJ[,1]))
annot <- annot[matchAnnot,]

# Throw away HJ probes where we dont have a annotation for
CTVHJ <- CTVHJ[which(as.character(CTVHJ[,1]) %in% as.character(annot[,6])),]
sortA <- match(CTVHJ[,1], annot[,6])

CTVHJ <- cbind(annot[sortA,5:6], CTVHJ)
# ONLY TAKE GOOD PROBES #

Affy <- read.csv("expHGU133A_CellType_Ratios.txt", row.names=1) # Load affy data
Affy[,2] <- as.character(Affy[,2])
Affy <- Affy[which(Affy[,2] != ''),] #Remove Non Annotated probes
cat("Affymetrix:", nrow(Affy), "Probes", length(unique(Affy[,2])), "genes\n")

Illu <- read.csv("expIllumina_CellType_Ratios.txt", row.names=1) # Load illumina data
Illu[,5] <- as.character(Illu[,5])
Illu <- removeMissing(Illu,5)
cat("Illumina:", nrow(Illu), "Probes", length(unique(Illu[,5])), "Unique genes\n")

Affy <- Affy[which(Affy[,2] %in% Illu[,5]),]                  # Match em up
Illu <- Illu[which(Illu[,5] %in% Affy[,2]),]                  # Match em up
CTVHJ <- CTVHJ[which(as.character(CTVHJ[,1]) %in% Illu[,5]),] # Match em up

cat("Affymetrix:", nrow(Affy), "Probes", length(unique(Affy[,2])), "Unique genes\n")
cat("Illumina:", nrow(Illu), "Probes", length(unique(Illu[,5])), "Unique genes\n")
cat("CTVHJ:", nrow(CTVHJ),"Probes", length(unique(CTVHJ[,1])), "Unique genes\n")

IlluVal <- 9:ncol(Illu) ; AffyVal <- 3:ncol(Affy) ; HJVal <- 4:ncol(CTVHJ) # Where is the data

colnames(Affy) <- paste0(colnames(Affy), "_A") # Add _A for Affy
colnames(Illu) <- paste0(colnames(Illu), "_I") # Add _I for Illumina

allGeneMatrix <- NULL                      # Create All Gene Vector Matrix index entrezGene
for(entrezGene in unique(Illu[,5])){
  probesI  <- which(Illu[,5]  == entrezGene) # Probe Indices in Illu
  probesA  <- which(Affy[,2]  == entrezGene) # Probe Indices in Affy
  probesHJ <- which(CTVHJ[,1] == entrezGene) # Probe Indices in CTVHJ

  IlluMeanProbe <- apply(t(t(Illu[probesI,IlluVal])), 2, function(x){mean(as.numeric(x))})
  if(any(is.na(IlluMeanProbe))) cat("gene:", entrezGene, "has NA in Illumina\n")
  
  AfftMeanProbe <- apply(t(t(Affy[probesA,AffyVal])), 2, function(x){mean(as.numeric(x))})
  if(any(is.na(AfftMeanProbe))) cat("gene:", entrezGene, "has NA in Affymetrics\n")

  HjMeanProbe <- apply(t(t(CTVHJ[probesHJ,HJVal])), 2, function(x){mean(as.numeric(x))})
  if(any(is.na(HjMeanProbe))) cat("gene:", entrezGene, "has NA in HJ gene\n")
  names(HjMeanProbe) <- c("Neurophil_HJ")

  row <- c(entrezGene, length(probesI), length(probesA), IlluMeanProbe, AfftMeanProbe, HjMeanProbe)
  allGeneMatrix <- rbind(allGeneMatrix, row)
}
rownames(allGeneMatrix) <- allGeneMatrix[,1]

write.csv(allGeneMatrix, file="cellT_Ratio_HUGO_Mean.txt", quote = FALSE)

# IF I take a more extreme change = less genes. Correlation goes up to 0.63
changed <- which(apply(allGeneMatrix[,7:ncol(allGeneMatrix)],1,function(x){any(abs(as.numeric(x)) > 10)}))

# Correlate every celltype against the other one
corM <- cor(apply(allGeneMatrix[,-c(1:3)],2,as.numeric), method="spearman", use="pair")
write.csv(corM, file="cellT_Cor_Spearman.txt", quote = FALSE)

cellT <- round(corM,3)[(1:length(IlluVal-2)),-c(1:length(IlluVal-2))]          # Print Illu Versus Correlations
cellT # Illumina vs Affy + Harm jan  

n <- colnames(corM) # After analysis to check the correlation of the corM against each other
for(x in (nrow(corM))){
  for(y in 1:(nrow(corM)-1)){
    pval <- cor.test(corM[,x],corM[,y])$p.value
    adjpval <- pval * 50; if(adjpval > 1){ adjpval = 1; }
    if(x != y){ #if(adjpval < 0.05){
      cat(n[x],"\t", n[y], "=", pval,"/", adjpval,"\n")
    }#}
  }
}

