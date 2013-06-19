# Downloads a GSE supplement dataset via FTP
downloadGSE <- function(GSE = "GSE14771"){
  tGSE <- paste0(strtrim(GSE, nchar(GSE)-3), "nnn")
  ftp  <- "ftp.ncbi.nlm.nih.gov"
  if(!file.exists(paste0(GSE,"_RAW.tar"))){
    system(paste0("wget ftp://",ftp,"/geo/series/",tGSE,"/", GSE,"/suppl/", GSE,"_RAW.tar"))
  }
  if(!file.exists(paste0(GSE,"_non-normalized.txt.gz"))){
    system(paste0("wget ftp://",ftp,"/geo/series/",tGSE,"/", GSE,"/suppl/",GSE,"_non-normalized.txt.gz"))
  }
  if(!file.exists(paste0(GSE,"_non-normalized_data.txt.gz"))){
    system(paste0("wget ftp://",ftp,"/geo/series/",tGSE,"/", GSE,"/suppl/",GSE,"_non-normalized_data.txt.gz"))
  }
  if(!file.exists(paste0(GSE))){ system(paste0("mkdir ", GSE)) }
  if(file.exists(paste0(GSE,"_RAW.tar"))){
    system(paste0("tar -xf ",GSE,"_RAW.tar -C ", GSE))
  }
  if(file.exists(paste0(GSE,"_non-normalized.txt.gz"))){
    system(paste0("gunzip -f ",GSE,"_non-normalized.txt.gz"))
    system(paste0("mv ",GSE,"_non-normalized.txt ", GSE))
  }
  if(file.exists(paste0(GSE,"_non-normalized_data.txt.gz"))){
    system(paste0("gunzip -f ",GSE,"_non-normalized_data.txt.gz"))
    system(paste0("mv ",GSE,"_non-normalized_data.txt ", GSE))
  }
}

# Download affymetrics aray data
downloadGSE("GSE22886") # 2 x 114 samples different cellTypes
downloadGSE("GSE6613")  # 105 samples, we use the 50 whole blood control
downloadGSE("GSE12288") # 222 samples
downloadGSE("GSE3846")  # 108
#downloadGSE("GSE1751")  # 31
#downloadGSE("GSE1343")  # 16
downloadGSE("GSE24250")  # 14

#AFFY GPL570
downloadGSE("GSE26440") # 130 samples
downloadGSE("GSE16059") # 88 samples

# Download illumina cellType specific array data
system(paste0("wget http://www.ebi.ac.uk/arrayexpress/files/E-TABM-633/E-TABM-633.processed.1.zip")
if(!file.exists("E-TABM-633")){ system(paste0("mkdir E-TABM-633")) } # Create folder
system(paste0("unzip E-TABM-633.processed.1.zip -d E-TABM-633")) # Extract file

# Download illumina array whole blood
downloadGSE("GSE24757")
downloadGSE("GSE19790")

