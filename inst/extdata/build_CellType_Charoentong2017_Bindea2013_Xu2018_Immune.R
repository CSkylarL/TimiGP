# Charoentong2017 Citation:https://www.sciencedirect.com/science/article/pii/S2211124716317090?via%3Dihub
# Charoentong, P., Finotello, F., Angelova, M., Mayer, C., Efremova, M., Rieder, D., ... & Trajanoski, Z. (2017). Pan-cancer immunogenomic analyses reveal genotype-immunophenotype relationships and predictors of response to checkpoint blockade. Cell reports, 18(1), 248-262.
# Download Table S6 and modify format: Charoentong2017_tableS6_Immune.txt


# Bindea2013 Citation:https://www.sciencedirect.com/science/article/pii/S1074761313004378?via%3Dihub
# Bindea, G., Mlecnik, B., Tosolini, M., Kirilovsky, A., Waldner, M., Obenauf, A. C., ... & Galon, J. (2013). Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity, 39(4), 782-795.
# Download Table S1 and modify format: Bindea2013_tableS1_immune.txt

# Xu2018 Citation:https://aacrjournals.org/cancerres/article/78/23/6575/543748/Xu2018-A-Web-Server-for-Resolving-Tumor
# Xu, L., Deng, C., Pang, B., Zhang, X., Liu, W., Liao, G., ... & Li, X. (2018). Xu2018: a web server for resolving tumor immunophenotype profiling. Cancer research, 78(23), 6575-6580.
# Download from http://biocc.hrbmu.edu.cn/Xu2018/download/signature%20annotation.txt: Xu2018_Signature.txt


## setwd("./") # working directory is the TimiGP folder
library(dplyr)
library(data.table)
rm(list = ls())
myinf1 <- "inst/extdata/Charoentong2017_tableS6_Immune.txt"
myinf2 <- "inst/extdata/Bindea2013_tableS1_immune.txt"
myinf3 <- "inst/extdata/Xu2018_Signature.txt"
myinf4 <- "inst/extdata/CellType_summary.csv"
myoutf1 <- "data/CellType_Charoentong2017_Bindea2013_Xu2018_Immune.rda"


## Charoentong2017 ----
myset1 <- read.table(myinf1, 
                     header = T, sep = "\t", stringsAsFactors = F)
head(myset1)
myset1 <- myset1[, -3]
colnames(myset1) <- c("Gene", "CellType")
myset1$Dataset <- "Charoentong2017"
myset1 <- myset1[, 3:1]
head(myset1)
sum(duplicated(myset1$Gene)) # 0

## Bindea2013 ----
myset2 <- read.table(myinf2, 
                     header = F, sep = "\t", stringsAsFactors = F)
head(myset2)
myset2 <- myset2[, -3]
colnames(myset2) <- c("CellType", "Gene")
myset2$Dataset <- "Bindea2013"
myset2 <- myset2[, c(3, 1,2)]
head(myset2)
se <- which(duplicated(myset2))
length(se)  #95
myset2 <- myset2[-se, ]
dim(myset2) # 586   3
se <- which(duplicated(myset2$Gene))
length(se)  # 7
myset2[se, ]
dup.gene <- myset2[se, 3]
dup.gene 
se <- which(myset2$Gene %in% dup.gene)
se # 38  44  45  55  65  74 115 231 276 283 314 334 579
myset2[se, ]
# don not remove duplicates(little influence)
# If keep mutually exclussive, you can keep the maker of subtype like below)------
# Dataset            CellType   Gene
# 44  Bindea2013             T cells   CD28 remove 38
# 53  Bindea2013             T cells    TRA remove 44
# 54  Bindea2013             T cells  ITM2A remove 45
# 67  Bindea2013      T helper cells  ITM2A
# 77  Bindea2013      T helper cells    TRA remove 65
# 87  Bindea2013      T helper cells   CD28
# 129 Bindea2013                 Tem    TRA 
# 264 Bindea2013         CD8 T cells  APBA2 remove 231
# 315 Bindea2013     Cytotoxic cells  APBA2 
# 323 Bindea2013            NK cells   XCL1 remove 283
# 356 Bindea2013            NK cells IGFBP5
# 381 Bindea2013 NK CD56bright cells   XCL1
# 674 Bindea2013       Normal mucosa IGFBP5 remove 579
# myset2 <- myset2[-se[c(1:3, 5, 8, 10, 13)], ]
# sum(duplicated(myset2$Gene)) #0
se <- which(myset2$CellType == "SW480 cancer cells")
length(se)
myset2_wo_cancer <- myset2[-se, ]
## Xu2018 ----

myset3 <- read.table(myinf3, 
                     header = T, sep = "\t", stringsAsFactors = F)
head(myset3)
myset3 <- myset3[, c(1,4)]
colnames(myset3) <- c("Gene", "CellType")
myset3$Dataset <- "Xu2018"
myset3 <- myset3[, 3:1]
head(myset3)

se <- which(myset3$CellType != "Multiple")
myset4 <- myset3[se,]
head(myset4)  # 0

se <- which(duplicated(myset4$Gene))
length(se)  # 38
length(se)/nrow(myset4)  # 0.4470588
table(myset4$CellType)
## Did not remove duplicated genes(too many)
## Some cell type without gene markers if remove them

## Combine annotation set(remove multiple in Xu2018) ---- 

geneset <- rbind(myset1, myset2_wo_cancer, myset4)
table(geneset$CellType,geneset$Dataset)
geneset$CellType <- paste0(geneset$Dataset,"_",geneset$CellType)
##Change cell type name
ref <- read.csv(myinf4,header = T)
refxx <- filter(ref,Dataset == "Charoentong2017"|Dataset == "Bindea2013"|Dataset == "Xu2018")
refxx$Paper.Name <- paste0(refxx$Dataset,"_",refxx$Paper.Name)
refxx <- refxx[,c("Method.Name", "Paper.Name")]
sum(!unique(geneset$CellType) %in% unique(refxx$Paper.Name)) #5
## unique(geneset$CellType)[!unique(geneset$CellType) %in% unique(refxx$Paper.Name)]
## pDC & Treg: 1 marker; Th 17: 3 markers; remove
## Only keep immune gene
## [1]  "Bindea2013_Th17 cells" "Bindea2013_TReg"          "Bindea2013_pDC"           "Bindea2013_Normal mucosa" "Bindea2013_Blood vessels" "Bindea2013_Lymph vessels"

sum(!unique(refxx$Paper.Name) %in% unique(geneset$CellType)) #1
# "Bindea2013_SW480 cancer cells"unique(refxx$Paper.Name)[!unique(refxx$Paper.Name) %in% unique(geneset$CellType)]
# "Bindea2013_SW480 cancer cells"
geneset <- merge(geneset,refxx,by.x = "CellType", by.y = "Paper.Name")

# "CellType" "Gene"     "Dataset" 
geneset <- geneset[4:2]
setnames(geneset,old = "Method.Name", new = "CellType")
dim(geneset) # 1375   3

CellType_Charoentong2017_Bindea2013_Xu2018_Immune <- geneset
# save(CellType_Charoentong2017_Bindea2013_Xu2018_Immune, file = myoutf1)


