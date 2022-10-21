# Tirosh2016 Citation:https://www.science.org/doi/full/10.1126/science.aad0501
# Tirosh, I., Izar, B., Prakadan, S. M., Wadsworth, M. H., Treacy, D., Trombetta, J. J., ... & Garraway, L. A. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science, 352(6282), 189-196.
# Download Table S3 and modify format: Tirosh2016_tableS3_melanoma.txt

## setwd("./") # working directory is the TimiGP folder
library(dplyr)
library(data.table)
rm(list = ls())
myinf1 <- "inst/extdata/Tirosh2016_tableS3_melanoma.txt"
myinf2 <- "inst/extdata/CellType_summary.csv"
myoutf1 <- "data/CellType_Tirosh2016_melanoma_TME.rda"


## Tirosh2016 #######
myset1 <- read.table(myinf1, 
                     header = T, sep = "\t", stringsAsFactors = F)
myset1$Dataset <- "Tirosh2016"

geneset <- myset1[c("Dataset", "CellType" ,"Gene")]


##Change cell type name
ref <- read.csv(myinf2,header = T)
refxx <- filter(ref,Dataset == "Tirosh2016")
refxx <- refxx[,c("Method.Name", "Paper.Name")]
sum(!unique(geneset$CellType) %in% unique(refxx$Paper.Name)) #0
sum(!unique(refxx$Paper.Name) %in% unique(geneset$CellType)) #0
geneset <- merge(geneset,refxx,by.x = "CellType", by.y = "Paper.Name")

geneset <- geneset[4:2]
setnames(geneset,old = "Method.Name", new = "CellType")
dim(geneset) # 663   3

table(geneset[1])

# B        CAF       Endo Macrophage   Melanoma          T 
# 31         88         95         92         47         38 

## save #######
CellType_Tirosh2016_melanoma_TME <- geneset
# save(CellType_Tirosh2016_melanoma_TME, file = myoutf1)