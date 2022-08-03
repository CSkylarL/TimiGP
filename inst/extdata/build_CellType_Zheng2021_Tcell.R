# Zheng2021 Citation:https://www.science.org/doi/10.1126/science.abe6474
# Zheng, L., Qin, S., Si, W., Wang, A., Xing, B., Gao, R., ... & Zhang, Z. (2021). Pan-cancer single-cell landscape of tumor-infiltrating T cells. Science, 374(6574), abe6474. 
# Download Table S2 and modify format: Zhang2021_tableS2_Tcell.txt


## setwd("./") # working directory is the TimiGP folder
rm(list = ls())
library(dplyr)

myinf1 <- "inst/extdata/Zheng2021_tableS2_Tcell.txt"
myoutf1 <- "data/CellType_Zheng2021_Tcell.rda"

myset1 <- read.table(myinf1,sep = "\t",header = T,stringsAsFactors = F)
myset1 <- myset1[,1:6]
colnames(myset1)
myset1$Gene <- paste(myset1[,3],myset1[,4],myset1[,5],myset1[,6] ,sep= ",")
head(myset1$Gene,1)
myset1 <- myset1[,c(1,2,7)]
myset1$cell <- myset1$meta.cluster %>% strsplit(split = ".c") %>% sapply("[[",1)

setS2 <- data.frame()

for (i in 1:nrow(myset1)) {
  xx <- strsplit(myset1$Gene[i],split = ",") %>% unlist() %>%
    sub(pattern = " ",replacement = "") 
  xx <- xx[xx!=""]
  nn <- length(xx)
  setS2 <- rbind(setS2,
                 data.frame(
                   Dataset=rep("Zheng2021",nn),
                   CellType = rep(paste0(myset1$cell[i],"+",myset1$Functional.properties[i]),nn),
                   Gene=xx,
                   Cluster =  rep(myset1$meta.cluster[i],nn),
                   Function = rep(myset1$Functional.properties[i],nn) ))
}

table(setS2$CellType)
length(unique(setS2$Gene)) # 389


# Many duplicates --- will not be removed
geneset <- setS2
# Remove duplicates
table(geneset$CellType)
se <- which(duplicated(geneset$Gene))
dup.gene <- geneset$Gene[se]
length(unique(dup.gene)) # 212

table(geneset$CellType)
colnames(geneset)
dim(geneset) #  926   5

geneset <- geneset[c(2,3,1)]
CellType_Zheng2021_Tcell <- geneset
# save(CellType_Zheng2021_Tcell, file = myoutf1)