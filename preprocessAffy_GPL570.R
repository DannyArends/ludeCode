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
GSE26440 <- read.csv("dataDescr/GSE26440.txt", sep="\t", header=FALSE, colClasses=("character")) # 130 samples

