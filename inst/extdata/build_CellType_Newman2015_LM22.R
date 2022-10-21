# Newman2015 Citation:https://www.nature.com/articles/nmeth.3337
# Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., ... & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature methods, 12(5), 453-457.# Download Table S6 and modify format: Newman2015_tableS6_Immune.txt
# Downloaded from table S1 and modified as: Newman2015_TableS1.txt

## setwd("./") # working directory is the TimiGP folder
library(dplyr)
library(data.table)
rm(list = ls())
myinf1 <- "inst/extdata/Newman2015_TableS1.txt"
myinf2 <- "inst/extdata/Aran2017_TableS3.txt"
myinf2.1 <- "inst/extdata/Aran2017_TableS1.csv"
myinf3 <- "inst/extdata/Luca2021_TableS4.txt"
myinf4 <- "inst/extdata/CellType_summary.csv"
myoutf1 <- "data/CellType_Newman2015_LM22.rda"


## Newman2015 ----
mydata <- read.table(myinf1, row.names = 1,
                     header = T, sep = "\t", stringsAsFactors = F)

myset1 <- data.frame()
for (i in 1:ncol(mydata)){
  se <- which(mydata[i] == 1)
  tmp <- data.frame(
                    CellType = colnames(mydata)[i],
                    Gene = rownames(mydata)[se],
                    Dataset = rep( "Newman2015",length(se)))
  myset1 <- rbind(myset1,tmp)
}




ref <- read.csv(myinf4,header = T)
refxx <- filter(ref,Dataset == "Newman2015")
refxx <- refxx[,c("Method.Name", "Paper.Name")]
sum(!unique(myset1$CellType) %in% unique(refxx$Paper.Name)) #0
sum(!unique(refxx$Paper.Name) %in% unique(myset1$CellType))#

myset1 <-  merge(myset1,refxx,by.x = "CellType", by.y = "Paper.Name")

myset1 <- myset1[c("Method.Name", "Gene",    "Dataset" )]
setnames(myset1,old = "Method.Name", new = "CellType")

CellType_Newman2015_LM22 <-  myset1
                                                 
#save(CellType_Newman2015_LM22, file = myoutf1)


