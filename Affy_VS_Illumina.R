#Match up

Illu <- read.csv("expIllumina_CellType_Ratios.txt", row.names=1)
Affy <- read.csv("expHGU133A_CellType_Ratios.txt", row.names=1)

Affy <- Affy[which(Affy[,2] %in% Illu[,2]),]
Illu <- Illu[which(Illu[,2] %in% Affy[,2]),]

cat("Affymetrix:", nrow(Affy), "Probes", length(unique(Affy[,2])), "Unique genes\n")
cat("Illumina:", nrow(Illu), "Probes", length(unique(Illu[,2])), "Unique genes\n")

IlluVal <- 3:ncol(Illu)
AffyVal <- 3:ncol(Affy)
colnames(Affy) <- paste0(colnames(Affy), "_A") # Add A for Affy
colnames(Illu) <- paste0(colnames(Illu), "_I") # Add I for Illumina

allGeneMatrix <- NULL                      # Create All Gene Vector Matrix index entrezGene
for(entrezGene in unique(Illu[,2])){
  probesI <- which(Illu[,2]==entrezGene)
  probesA <- which(Affy[,2]==entrezGene)
  IlluMeanProbe <- apply(t(t(Illu[probesI,IlluVal])), 2, mean)
  AfftMeanProbe <- apply(t(t(Affy[probesA,AffyVal])), 2, mean)
  row <- c(entrezGene, length(probesI), length(probesA), IlluMeanProbe, AfftMeanProbe)
  allGeneMatrix <- rbind(allGeneMatrix, row)
}

write.csv(allGeneMatrix, file="cellT_Ratio_Per_Gene.txt", quote = FALSE)

# Correlate every celltype against the other ones
corM <- cor(allGeneMatrix[,-c(1:3)], method="spearman", use="pair")
cellT <- round(corM,3)[(IlluVal-2),-c(IlluVal-2)]          # Print Illu Versus Correlations
write.csv(cellT, file="cellT_Cor_Spearman.txt", quote = FALSE)
