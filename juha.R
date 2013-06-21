wholeblood <- read.csv("GPL570_WholeBlood.txt",sep='\t',row.names=1)
neutrophil <- read.csv("GPL570_Neutrophil.txt",sep='\t',row.names=1)

badneutros <- c("GSM141250", "GSM141251", "GSM141252", "GSM141253", "GSM141254", "GSM141255", "GSM141256",
                "GSM141257", "GSM549581", "GSM549582", "GSM549583", "GSM549584")

neutrophil <- neutrophil[,-which(colnames(neutrophil) %in% badneutros)]

zs <- NULL
for(x in 1:nrow(wholeblood)){
  p  <- t.test(wholeblood[x,],neutrophil[x,])$p.value
  zs <- c(zs, qnorm(1 - (p/2)))
}

