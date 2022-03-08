# Zhang2021S2 Citation:https://www.science.org/doi/10.1126/science.abe6474
# Zheng, L., Qin, S., Si, W., Wang, A., Xing, B., Gao, R., ... & Zhang, Z. (2021). Pan-cancer single-cell landscape of tumor-infiltrating T cells. Science, 374(6574), abe6474. 
# Download Table S2 and modify format: Zhang2021_tableS2_Tcell.txt


## setwd("./") # working directory is the TimiGP folder
rm(list = ls())
library(dplyr)

myinf1 <- "inst/extdata/Zhang2021_tableS2_Tcell.txt"
myoutf1 <- "data/CellType_Zhang2021S2_Tcell.rda"

myset1 <- read.table(myinf1,sep = "\t",header = T,stringsAsFactors = F)
myset1 <- myset1[,1:6]
colnames(myset1)
myset1$Gene <- paste(myset1[,3],myset1[,4],myset1[,5],myset1[,6] ,sep= ",")
head(myset1$Gene,1)
myset1 <- myset1[,c(1,2,7)]

setS2 <- data.frame()

for (i in 1:nrow(myset1)) {
  xx <- strsplit(myset1$Gene[i],split = ",") %>% unlist() %>%
    sub(pattern = " ",replacement = "") 
  xx <- xx[xx!=""]
  nn <- length(xx)
  setS2 <- rbind(setS2,
                 data.frame(
                   Dataset=rep("Zhang2021S2",nn),
                   CellType = rep(paste0(myset1$meta.cluster[i],"(",myset1$Functional.properties[i],")"),nn),
                   Gene=xx,
                   Cluster =  rep(myset1$meta.cluster[i],nn),
                   Function = rep(myset1$Functional.properties[i],nn) ))
}

table(setS2$CellType)
length(unique(setS2$Gene)) # 389


# Too many duplicates
geneset <- setS2
# Remove duplicates
table(geneset$CellType)
se <- which(duplicated(geneset$Gene))
dup.gene <- geneset$Gene[se]
length(unique(dup.gene)) # 212
# dup.gene 
# se <- which(geneset$Gene %in% dup.gene)
# length(se)
# geneset <- geneset[-se,]
table(geneset$CellType)
colnames(geneset)
dim(geneset) #  926   5

geneset <- geneset[c(2,3,1)]
CellType_Zhang2021S2_Tcell <- geneset
# save(CellType_Zhang2021S2_Tcell, file = myoutf1)