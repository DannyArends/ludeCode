setwd("/home/danny/Github/LudeNew")

ProbeAnnotation <- read.csv("GPL6102-11574.txt",sep="\t",skip=27) # Probe annotation

CellTypeDATA <- read.table("E-TABM-633/HaemAtlasMKEBNormalizedIntensities.csv",sep='\t',row.names=1,header=TRUE)
CellTypeDATA <- CellTypeDATA[-1,]
tmp <- apply(CellTypeDATA, 2, as.numeric)
rownames(tmp) <- rownames(CellTypeDATA)
CellTypeDATA  <- tmp

CellTypeAnnotation <- read.csv("dataDescr/E-TABM-633.txt",sep="\t", header=TRUE)
cellTypeNames <- as.character(CellTypeAnnotation[, "Hybridization.Name"])
cellTypeTypes <- as.character(CellTypeAnnotation[, "Extract.Name"])
cellTypes     <- unlist(lapply(strsplit(cellTypeTypes,"-"),"[",3))

cellTypeAnnot <- cbind(cellTypeNames, cellTypes) # FInally the annotation we want Hyb ref -> Cell-type
CellTypeDATA  <- CellTypeDATA[, match(cellTypeAnnot[,1], gsub("X","",colnames(CellTypeDATA)))] # ARRANGE

colnames(CellTypeDATA) <- as.character(cellTypeAnnot[,2])

ProbeAnnotation <- ProbeAnnotation[which(ProbeAnnotation[,"Symbol"] != ""),]

IlluMean <- apply(CellTypeDATA[,which(colnames(CellTypeDATA) == "Granulocyte")],1,mean)

inAnnot <- which(names(IlluMean) %in% ProbeAnnotation[,1])
IlluMean <- IlluMean[inAnnot]
sortAnnot <- match(names(IlluMean), ProbeAnnotation[,1]) # Align
MeanMatrix <- cbind(as.character(ProbeAnnotation[sortAnnot,"Symbol"]), IlluMean)

setwd("~/Github/Juha/")
Affy <- read.csv("Neutrophil_RNAseq.txt",sep='\t')

inMean <- which(MeanMatrix[,1] %in% rownames(Affy))
inAffy <- which(rownames(Affy) %in% MeanMatrix[,1])
MeanMatrix <- MeanMatrix[inMean, ]
Affy <- Affy[inAffy, ]

RnaAffyIllu <-NULL
for(gene in rownames(Affy)){
  illumean <- mean(as.numeric(MeanMatrix[which(MeanMatrix[,1] == gene),2]))
  AffyRowID <- which(rownames(Affy)==gene)
  RnaAffyIllu <- rbind(RnaAffyIllu, c(gene, illumean, Affy[AffyRowID,-8]))
}




