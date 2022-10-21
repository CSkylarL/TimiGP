# Newman2015 Citation:https://www.nature.com/articles/nmeth.3337
# Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., ... & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature methods, 12(5), 453-457.# Download Table S6 and modify format: Newman2015_tableS6_Immune.txt
# Downloaded from table S1 and modified as: Newman2015_TableS1.txt

# Aran2017 Citation:https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1
# Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome biology, 18(1), 1-14.
# Downloaded from table S3 and modified as: Aran2017_TableS3.txt


#Luca2021  Citation:https://www.sciencedirect.com/science/article/pii/S0092867421010618?via%3Dihub
# Luca, B. A., Steen, C. B., Matusiak, M., Azizi, A., Varma, S., Zhu, C., ... & Newman, A. M. (2021). Atlas of clinically distinct cell states and ecosystems across human solid tumors. Cell, 184(21), 5482-5496.# Download Table S1 and modify format: Aran2017_tableS1_immune.txt
# Downloadedfrom table S4 and modified as: Luca2021_TableS4.txt

## setwd("./") # working directory is the TimiGP folder
library(dplyr)
library(data.table)
rm(list = ls())
myinf1 <- "inst/extdata/Newman2015_TableS1.txt"
myinf2 <- "inst/extdata/Aran2017_TableS3.txt"
myinf2.1 <- "inst/extdata/Aran2017_TableS1.csv"
myinf3 <- "inst/extdata/Luca2021_TableS4.txt"
myinf4 <- "inst/extdata/CellType_summary.csv"
myoutf1 <- "data/CellType_Newman2015_Aran2017_Luca2021_list.rda"


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

## Aran2017 ----
mydata <- read.table(myinf2, 
                     header = T, sep = "\t", stringsAsFactors = F)

mydata$celltype <- mydata$Celltype_Source_ID %>%
  strsplit("_",fixed = T) %>%
  lapply("[[", 1) %>%
  unlist()

cid <- mydata$celltype %>% unique()


myset2 <- data.frame()
for (i in cid){
  se <- which(mydata$celltype == i)
  mygen <- mydata[se,]$Gene %>% strsplit("|",fixed = T) %>% unlist() %>% unique()
  tmp <- data.frame(
    CellType = rep(i,length(mygen)),
    Gene = mygen,
    Dataset = rep( "Aran2017",length(mygen)))
  myset2 <- rbind(myset2,tmp)
}

anno <- read.csv(myinf2.1)

sum(!unique(myset2$CellType) %in% unique(anno$Cell.types)) #0
sum(!unique(anno$Cell.types) %in% unique(myset2$CellType))#0

myset2 <- merge(anno,myset2,by.x="Cell.types",by.y="CellType")


ref <- read.csv(myinf4,header = T)
refxx <- filter(ref,Dataset == "Aran2017")
refxx <- refxx[,c("Method.Name", "Paper.Name")]
sum(!unique(myset2$Cell.types) %in% unique(refxx$Paper.Name)) #0
sum(!unique(refxx$Paper.Name) %in% unique(myset2$Cell.types))#0
myset2 <- merge(refxx,myset2,by.x="Paper.Name",by.y="Cell.types")

myset2 <- myset2 [c("Method.Name" , "Gene"  ,
                      "Dataset" ,"Group"     , "Subgroup"       )]

setnames(myset2,old = "Method.Name", new = "CellType")

# Luca2021 -----

myset3 <- read.table(myinf3, 
                     header = T, sep = "\t", stringsAsFactors = F)




ref <- read.csv(myinf4,header = T)
refxx <- filter(ref,Dataset == "Luca2021")
refxx <- refxx[,c("Method.Name", "Paper.Name")]
sum(!unique(myset3$Cell.type) %in% unique(refxx$Paper.Name)) #0
sum(!unique(refxx$Paper.Name) %in% unique(myset3$Cell.type))#

myset3 <-  merge(myset3,refxx,by.x = "Cell.type", by.y = "Paper.Name")

myset3$Dataset <- "Luca2021"
myset3$CellType <- paste0(myset3$Method.Name,".",myset3$Cell.state)

myset3 <- myset3[c("CellType", "Gene",    "Dataset" ,  "Method.Name","Cell.state")]
setnames(myset3,old = "Method.Name", new = "Cell")
setnames(myset3,old = "Cell.state" , new = "State")


CellType_Newman2015_Aran2017_Luca2021_list <- list(Newman2015 = myset1,
                                                   Aran2017 = myset2,
                                                   Luca2021 = myset3)
#save(CellType_Newman2015_Aran2017_Luca2021_list, file = myoutf1)


