source("http://bioconductor.org/biocLite.R")    # Connect to bioConductor
#biocLite(c("affy","lumi"))                      # Install the affy package
# ILLUMINA #
# cellTypes Samples: HaemAtlasMKEBNormalizedIntensities.csv
# wholeBlood Samples: GSE24757 and GSE19790
library(lumi)
library(preprocessCore)

readIlluminaDir <- function(folder, annotation){
  data <- NULL
  for(x in paste0(folder, "/", annotation[,1], ".txt")){
    cat(x,"\n")
    m <- read.csv(x, header=TRUE, row.names=1,sep="\t")
    data <- cbind(data, apply(m,2,as.numeric)[,2])
  }
  rownames(data) <- rownames(m)
  data <- cbind(m[,1], data)
  colnames(data) <- c("Array_Address_Id", annotation[,1])
  invisible(data)
}

ProbeAnnotation <- read.csv("GPL6102-11574.txt",sep="\t",skip=27) # Probe annotation
GSE24757A <- read.csv("dataDescr/GSE24757.txt", sep="\t", header=FALSE, colClasses=("character")) # Annotation
GSE19790A <- read.csv("dataDescr/GSE19790.txt", sep="\t", header=FALSE, colClasses=("character")) # Annotation
GSE13255A <- read.csv("dataDescr/GSE13255.txt", sep="\t", header=FALSE, colClasses=("character")) # Annotation

GSE19790DATA    <- read.csv("GSE19790/GSE19790_non-normalized.txt",sep="\t",row.names=1)
GSE19790DATA    <- GSE19790DATA[, grep("AVG_Signal", colnames(GSE19790DATA))]
colnames(GSE19790DATA) <- gsub(".AVG_Signal","", colnames(GSE19790DATA))
colnames(GSE19790DATA) <- gsub("X","Ind", colnames(GSE19790DATA))

GSE13255DATA <- readIlluminaDir("GSE13255", GSE13255A)

GSE24757DATA <- read.csv("GSE24757/GSE24757_non-normalized_data.txt",sep="\t",row.names=1, skip=4,header=TRUE)

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

GSE24757DATA <- GSE24757DATA[match(rownames(CellTypeDATA), rownames(GSE24757DATA)), ] # Align data to CellType
GSE19790DATA <- GSE19790DATA[match(rownames(CellTypeDATA), rownames(GSE19790DATA)), ] # Align data to CellType
GSE13255DATA <- GSE13255DATA[which(rownames(GSE13255DATA) %in% rownames(CellTypeDATA)),]
GSE13255DATA <- GSE13255DATA[match(rownames(CellTypeDATA), rownames(GSE13255DATA)),]  # Align data to CellType

CellTypeDATA[1:5, 1:5]
CellTypeAnnotation[1:5, 1:5]
GSE19790DATA[1:5, 1:5]
GSE24757DATA[1:5, 1:5]
GSE13255DATA[1:5, 1:5]

phenotypes <- unique(cellTypeAnnot[,2])  # Calculate CellType Means
cellMeans <- NULL
for(phe in phenotypes){
  ctdIndices <- which(cellTypeAnnot[,2] == phe)
  cellMeans <- cbind(cellMeans, apply(CellTypeDATA[, ctdIndices], 1, mean))
}
rownames(cellMeans) <- rownames(CellTypeDATA)
colnames(cellMeans) <- phenotypes

if(any(!(rownames(GSE19790DATA) == rownames(GSE24757DATA)))){ # Check if the rownames MATCH
  stop("!!!!!")
}
wholeBlood <- cbind(GSE19790DATA, GSE24757DATA, GSE13255DATA[, -1])
#write.csv(wholeBlood, file="expIllumina_WholeBlood_RAW.txt", quote = FALSE)

#QUANTILE Normalization and log2 transform makes it equal to cellMeans
wholeBlood <- normalize.quantiles(apply(wholeBlood, 2, as.numeric))
wholeBlood <- apply(wholeBlood,2, log2)
rownames(wholeBlood) <- rownames(GSE19790DATA)
colnames(wholeBlood) <- c(colnames(GSE19790DATA), colnames(GSE24757DATA), colnames(GSE13255DATA)[-1])
#write.csv(wholeBlood, file="expIllumina_WholeBlood_QNORM_Log2.txt", quote = FALSE)

wholeMean <- apply(wholeBlood, 1, mean)
wholeSd <- apply(wholeBlood, 1, sd, na.rm=TRUE)

cellMeans <- cellMeans[which(wholeSd!=0),] #Remove SD==0
wholeMean <- wholeMean[which(wholeSd!=0)]
wholeSd <- wholeSd[which(wholeSd!=0)]

cellTypeRatios <- (cellMeans - wholeMean) / wholeSd #cellMeans/wholeMean

hasID <- which(!is.na(ProbeAnnotation[,15]))
ProbeAnnotation <- ProbeAnnotation[hasID,]

sortG <- match(rownames(cellTypeRatios), ProbeAnnotation[,1])
annotatedRatios <- cbind(ProbeAnnotation[sortG,c(1,15)],cellTypeRatios)

annot <- read.table("2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable-H8v2AndHT12v3.txt",sep='\t',header=TRUE)

CellTypeVectorHJ <- read.csv("CellTypeSpecificityMatrix.txt", sep='\t')
CellTypeVectorHJ <- CellTypeVectorHJ[-c(8229,8230),]
annot <- annot[which(annot[,6] %in% CellTypeVectorHJ[,1]),]

annot <- annot[which(annot[,15] %in% annotatedRatios[,2]), c(1:5,15)] # Take only our annotations
matchA <- which(annotatedRatios[,2] %in% annot[,6])
annotatedRatios <- annotatedRatios[matchA,]

sortA <- match(annotatedRatios[,2], annot[,6])
Hugo <- cbind(annot[sortA,], annotatedRatios)
write.csv(cbind(annot[sortA,], annotatedRatios), file="expIllumina_CellType_Ratios.txt", quote = TRUE)

