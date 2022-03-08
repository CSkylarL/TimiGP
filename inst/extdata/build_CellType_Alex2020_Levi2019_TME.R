# Alex2020 Citation: https://www.sciencedirect.com/science/article/pii/S107476132030409X#app2
# Hughes, T. K., Wadsworth II, M. H., Gierahn, T. M., Do, T., Weiss, D., Andrade, P. R., ... & Shalek, A. K. (2020). Second-strand synthesis-based massively parallel scRNA-seq reveals cellular states and molecular features of human inflammatory skin pathologies. Immunity, 53(4), 878-894.
# Download Table S3 and modify format: Alex2020_tableS3_skin.txt


# Levi2019 Citation:https://www.science.org/doi/full/10.1126/science.aad0501
# Tirosh, I., Izar, B., Prakadan, S. M., Wadsworth, M. H., Treacy, D., Trombetta, J. J., ... & Garraway, L. A. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science, 352(6282), 189-196.
# Download Table S3 and modify format: Levi2019_tableS3_melanoma.txt

## setwd("./") # working directory is the TimiGP folder
library(dplyr)
rm(list = ls())
myinf1 <- "inst/extdata/Alex2020_tableS3_skin.txt"
myinf2 <- "inst/extdata/Levi2019_tableS3_melanoma.txt"
myoutf1 <- "data/CellType_Alex2020_Levi2019_TME.rda"


## Alex2020 ######
  myset1 <- read.table(myinf1, 
                       header = T, sep = "\t", stringsAsFactors = F)
  myset1$Dataset <- "Alex2020"
  
  # Seurat FindMarkers
  # p_val : p_val (unadjusted)
  # avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
  # pct.1 : The percentage of cells where the feature is detected in the first group
  # pct.2 : The percentage of cells where the feature is detected in the second group
  # p_val_adj : Adjusted p-value, based on bonferroni correction using all features in the dataset.
  myset1$min.diff.pct <- myset1$pct.1-myset1$pct.2
  
  unique(myset1$Gene) %>% length() # 7817
  table(myset1$CellType)
  colnames(myset1)
  
  # select markers
  myset1 <-myset1 %>%
    filter(p_val_adj<0.001 & avg_logFC >0 & min.diff.pct >0) %>%
    group_by(CellType) %>%
    arrange(-avg_logFC,-min.diff.pct,.by_group = TRUE) %>%
    slice(1:30) %>%
    ungroup() %>%
    data.frame()
  table(myset1$CellType)
  
  
  # Remove duplicates
  table(myset1$CellType)
  se <- which(duplicated(myset1$Gene))
  dup.gene <- myset1$Gene[se]
  dup.gene 
  se <- which(myset1$Gene %in% dup.gene)
  length(se)
  myset1 <- myset1[-se,]
  table(myset1$CellType)
  colnames(myset1)
  dim(myset1) #367   9
  
  
  
  # Dataset CellType Gene
  colnames(myset1)
  myset1 <- myset1[c("Dataset", "CellType" ,"Gene")]
  geneset <- myset1






## Levi2019 #######
myset1 <- read.table(myinf2, 
                     header = T, sep = "\t", stringsAsFactors = F)
myset1$Dataset <- "Levi2019"

myset1 <- myset1[c("Dataset", "CellType" ,"Gene")]


# add CAFs and melanoma  of Levi 2019 to Alex2020
# CAFs have many genes the same as Fibro(Alex2020)
# melanoma have many genes the same as melanocyt(Alex2020)
# macrophage --> TIL-m
# T cell --> TIL-T
# B cell --> TIL-B
myset1$CellType[which(myset1$CellType == "B cells")] <- "TIL B cells"
myset1$CellType[which(myset1$CellType == "T cells")] <- "TIL T cells"
myset1$CellType[which(myset1$CellType == "Macrophages")] <- "TIL Macrophages"

## combine two sets####

se <- which(myset1$CellType %in% c("CAFs","melanoma","TIL Macrophages",
                                   "TIL B cells","TIL T cells"))
geneset <- rbind(myset1[se,],geneset)

table(geneset$CellType)
# B            CAFs            Endo           Fibro    HairFollicle              KC      Langerhans       Lymphatic 
# 21              88              26              30              22              22              19              25 
# Mast      Melanocyte        melanoma         Myeloid          Plasma         Schwann        Sebocyte               T 
# 27              26              47              18              20              27              25              29 
# TIL B cells TIL Macrophages     TIL T cells            VSMC 
# 31              92              38              30 

dim(geneset) # 663  3

## format ######

geneset <- geneset[c(2,3,1)]
se <- which(geneset$CellType == "B")
geneset[se,1] <- "Skin B cells"
se <- which(geneset$CellType == "T")
geneset[se,1] <- "Skin T cells"
se <- which(geneset$CellType == "Plasma")
geneset[se,1] <- "Skin Plasma cells"
se <- which(geneset$CellType == "Myeloid")
geneset[se,1] <- "Skin Myeloid cells"
se <- which(geneset$CellType == "Mast")
geneset[se,1] <- "Skin Mast cells"
se <- which(geneset$CellType == "Langerhans")
geneset[se,1] <- "Langerhans cells"
se <- which(geneset$CellType == "Endo")
geneset[se,1] <- "Venular endothelial cells"
se <- which(geneset$CellType == "Fibro")
geneset[se,1] <- "Fibroblasts"
se <- which(geneset$CellType == "HairFollicle")
geneset[se,1] <- "Hair follicles"
se <- which(geneset$CellType == "KC")
geneset[se,1] <- "Keratinocytes"
se <- which(geneset$CellType == "Lymphatic")
geneset[se,1] <- "Lymphatic endothelial cells"
se <- which(geneset$CellType == "Melanocyte")
geneset[se,1] <- "Melanocytes"
se <- which(geneset$CellType == "Schwann")
geneset[se,1] <- "Schwann cells"
se <- which(geneset$CellType == "Sebocyte")
geneset[se,1] <- "Sebocytes"
se <- which(geneset$CellType == "VSMC")
geneset[se,1] <- "Vascular smooth muscle cells"
se <- which(geneset$CellType == "melanoma")
geneset[se,1] <- "Melanoma"
se <- which(geneset$CellType == "CAFs")
geneset[se,1] <- "Cancer-associated fibroblasts"
table(geneset[1])
## save #######
CellType_Alex2020_Levi2019_TME <- geneset
# save(CellType_Alex2020_Levi2019_TME, file = myoutf1)