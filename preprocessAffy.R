source("http://bioconductor.org/biocLite.R")   # Connect to bioConductor
#install.packages(c("preprocessCore", "biomaRt")
#biocLite(c("affy", "gcrma", "lumi", "gcrma","AnnotationDbi")) # Install the affy package
library(affy)                                                  # Load the affy package
library(gcrma)                                                 # gcRMA normalisation
library(lumi)                                                  # Load the illumina package
library(biomaRt)                                               # Biomart annotations
library(AnnotationDbi)                                         # Load the annotations package
library(preprocessCore)                                        # Normalization

# AFFYMETRIX #
# We need the array info files: (CELname \t Other Columns)
GSE22886A <- read.csv("dataDescr/GSE22886.txt", sep="\t", header=FALSE, colClasses=("character"))
GSE6613A <- read.csv("dataDescr/GSE6613.txt", sep="\t", header=FALSE, colClasses=("character"))

getCol <- function(x, aType){ GSE22886A[which(GSE22886A[,4]==aType), x] } # Helper function
controlSamples  <- GSE6613A[which(grepl("control", GSE6613A[,2])),1]# Only control samples, Why addcomplexity

HGU133A <- paste0("GSE22886/",  getCol(1,"[HG-U133A]"),".CEL.gz") # The Filenames for U133A
HGU133Awb <- paste0("GSE6613/", controlSamples,".CEL.gz")# Filenames for whole blood controls

# Loads CEL.gz data files and return nolog gcRMA expression data
loadCELdata <- function(filenames, file="", doQuantiles = FALSE){
  library(affy)
  data <- ReadAffy(filenames = filenames)
  eset <- gcrma(data)
  expr <- exprs(eset)
  arrayNames <- colnames(expr)
  probeNames <- rownames(expr)
  if(doQuantiles){
    expr <- normalize.quantiles(expr)
    colnames(expr) <- arrayNames
    rownames(expr) <- probeNames
  }
  if(file != "") write.csv(expr, file=file, quote = FALSE)
  invisible(expr)
}

# Calculate the cell type means for array type aType
calcCellTypeMeans <- function(expData, aType, file = ""){
  phenames   <- unique(getCol(2, aType))
  cTypeMeans <- NULL
  for(phe in phenames){
    cTypeMeans <- cbind(cTypeMeans, apply(expData[,which(getCol(2,aType) == phe)], 1, mean))
  }
  colnames(cTypeMeans) <- phenames
  if(file != "") write.csv(cTypeMeans, file=file, quote = FALSE)
  invisible(cTypeMeans)
}

if(!file.exists("expHGU133A_gcRMA_ALL.txt")){ # Load .CEL.gz data for the Affy HGU133A samples
  expression <- loadCELdata(c(HGU133A, HGU133Awb), "expHGU133A_gcRMA_ALL.txt", doQuantiles = TRUE)
}
expression <- read.csv("expHGU133A_gcRMA_ALL.txt", row.names=1)

wholeBlood <- expression[,gsub("GSE6613/", "", HGU133Awb)]
cellTypes  <- expression[,gsub("GSE22886/","", HGU133A)]

write.csv(wholeBlood, file="expHGU133A_gcRMA_WholeBlood.txt", quote = FALSE)
write.csv(cellTypes, file="expHGU133A_gcRMA_CellType.txt", quote = FALSE)

wholeMean <- apply(wholeBlood, 1, mean)
wholeSd <- apply(wholeBlood, 1, sd, na.rm=TRUE)

cellMeans <- calcCellTypeMeans(cellTypes, "[HG-U133A]", "expHGU133A_gcRMA_CellType_mean.txt")

cellMeans <- cellMeans[which(wholeSd!=0),] #Remove SD==0
wholeMean <- wholeMean[which(wholeSd!=0)]
wholeSd <- wholeSd[which(wholeSd!=0)]

#cellTypeRatios <- cellTypeRatios[which(wholeMean > 4),]
#wholeMean <- wholeMean[which(wholeMean > 4)]

cellTypeRatios <- (cellMeans - wholeMean) / wholeSd #cellMeans/wholeMean
write.csv(cellTypeRatios, file="expHGU133A_CellType_Ratios.txt", quote = FALSE)

ids    <- rownames(cellTypeRatios)
mart   <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
toGene <- getBM(attributes = c("affy_hg_u133a", "hgnc_symbol"),
                filters="affy_hg_u133a", values=ids, mart=mart)

toGene <- toGene[which(!is.na(toGene[,2])),] # Removed the ones in togene which don't have a gene
toGene <- toGene[which(!toGene[,2]==""),] # Removed the ones in togene which don't have a gene

matched <- which(rownames(cellTypeRatios) %in% toGene[,1])
cellTypeRatios <- cellTypeRatios[matched, ] # Removed the ones in cellTypeRatios which don't have a match

sortG <- match(rownames(cellTypeRatios), toGene[,1]) # Sorter to put toGene in the order of cellTypeRatios

write.csv(cbind(toGene[sortG,],cellTypeRatios), file="expHGU133A_CellType_Ratios.txt", quote = FALSE)

