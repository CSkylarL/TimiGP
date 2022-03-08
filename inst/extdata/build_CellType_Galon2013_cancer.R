# Galon2013 Citation:https://www.sciencedirect.com/science/article/pii/S1074761313004378?via%3Dihub
# Bindea, G., Mlecnik, B., Tosolini, M., Kirilovsky, A., Waldner, M., Obenauf, A. C., ... & Galon, J. (2013). Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity, 39(4), 782-795.
# Download Table S1 and modify format: Galon2013_tableS1_immune.txt


## setwd("./") # working directory is the TimiGP folder
rm(list = ls())
myinf2 <- "inst/extdata/Galon2013_tableS1_immune.txt"
myoutf1 <- "data/CellType_Galon2013_cancer.rda"


## Galon2013 ----
myset2 <- read.table(myinf2, 
                     header = F, sep = "\t", stringsAsFactors = F)
head(myset2)
myset2 <- myset2[, -3]
colnames(myset2) <- c("CellType", "Gene")
myset2$Dataset <- "Galon2013"
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
# Dataset            CellType   Gene
# 44  Galon2013             T cells   CD28 remove 38
# 53  Galon2013             T cells    TRA remove 44
# 54  Galon2013             T cells  ITM2A remove 45
# 67  Galon2013      T helper cells  ITM2A
# 77  Galon2013      T helper cells    TRA remove 65
# 87  Galon2013      T helper cells   CD28
# 129 Galon2013                 Tem    TRA 
# 264 Galon2013         CD8 T cells  APBA2 remove 231
# 315 Galon2013     Cytotoxic cells  APBA2 
# 323 Galon2013            NK cells   XCL1 remove 283
# 356 Galon2013            NK cells IGFBP5
# 381 Galon2013 NK CD56bright cells   XCL1
# 674 Galon2013       Normal mucosa IGFBP5 remove 579
myset2 <- myset2[-se[c(1:3, 5, 8, 10, 13)], ]
sum(duplicated(myset2$Gene)) #0


geneset <- myset2
geneset <- geneset[c(2,3,1)]

CellType_Galon2013_cancer <- geneset
# save(CellType_Galon2013_cancer, file = myoutf1)
