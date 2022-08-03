# Bindea2013 Citation:https://www.sciencedirect.com/science/article/pii/S1074761313004378?via%3Dihub
# Bindea, G., Mlecnik, B., Tosolini, M., Kirilovsky, A., Waldner, M., Obenauf, A. C., ... & Galon, J. (2013). Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity, 39(4), 782-795.
# Download Table S1 and modify format: Bindea2013_tableS1_immune.txt


## setwd("./") # working directory is the TimiGP folder
rm(list = ls())
myinf2 <- "inst/extdata/Bindea2013_tableS1_immune.txt"
myinf4 <- "inst/extdata/CellType_summary.csv"
myoutf1 <- "data/CellType_Bindea2013_cancer.rda"


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


geneset <- myset2



ref <- read.csv(myinf4,header = T)
refxx <- filter(ref,Dataset == "Bindea2013")
refxx <- refxx[,c("Method.Name", "Paper.Name")]
sum(!unique(geneset$CellType) %in% unique(refxx$Paper.Name)) #5
## unique(geneset$CellType)[!unique(geneset$CellType) %in% unique(refxx$Paper.Name)]
## pDC & Treg: 1 marker; Th 17: 3 markers; remove
## Only keep immune gene
## [1]  "Bindea2013_Th17 cells" "Bindea2013_TReg"          "Bindea2013_pDC"           "Bindea2013_Normal mucosa" "Bindea2013_Blood vessels" "Bindea2013_Lymph vessels"

sum(!unique(refxx$Paper.Name) %in% unique(geneset$CellType)) #0

geneset <- merge(geneset,refxx,by.x = "CellType", by.y = "Paper.Name")


# "CellType" "Gene"     "Dataset" 
geneset <- geneset[4:2]
setnames(geneset,old = "Method.Name", new = "CellType")
dim(geneset) # 545   3


CellType_Bindea2013_cancer <- geneset
# save(CellType_Bindea2013_cancer, file = myoutf1)
