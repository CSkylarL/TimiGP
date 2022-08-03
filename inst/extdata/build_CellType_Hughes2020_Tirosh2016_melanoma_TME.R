# Hughes2020 Citation: https://www.sciencedirect.com/science/article/pii/S107476132030409X#app2
# Hughes, T. K., Wadsworth II, M. H., Gierahn, T. M., Do, T., Weiss, D., Andrade, P. R., ... & Shalek, A. K. (2020). Second-strand synthesis-based massively parallel scRNA-seq reveals cellular states and molecular features of human inflammatory skin pathologies. Immunity, 53(4), 878-894.
# Download Table S3 and modify format: Hughes2020_tableS3_skin.txt


# Tirosh2016 Citation:https://www.science.org/doi/full/10.1126/science.aad0501
# Tirosh, I., Izar, B., Prakadan, S. M., Wadsworth, M. H., Treacy, D., Trombetta, J. J., ... & Garraway, L. A. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science, 352(6282), 189-196.
# Download Table S3 and modify format: Tirosh2016_tableS3_melanoma.txt

## setwd("./") # working directory is the TimiGP folder
library(dplyr)
library(data.table)
rm(list = ls())
myinf1 <- "inst/extdata/Hughes2020_tableS3_skin.txt"
myinf2 <- "inst/extdata/Tirosh2016_tableS3_melanoma.txt"
myinf3 <- "inst/extdata/CellType_summary.csv"
myoutf1 <- "data/CellType_Hughes2020_Tirosh2016_melanoma_TME.rda"


## Hughes2020 ######
  myset1 <- read.table(myinf1, 
                       header = T, sep = "\t", stringsAsFactors = F)
  myset1$Dataset <- "Hughes2020"
  
  # Seurat FindMarkers
  # p_val : p_val (unadjusted)s
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


## Tirosh2016 #######
myset1 <- read.table(myinf2, 
                     header = T, sep = "\t", stringsAsFactors = F)
myset1$Dataset <- "Tirosh2016"

myset1 <- myset1[c("Dataset", "CellType" ,"Gene")]


# add CAFs and melanoma  of Levi 2019 to Hughes2020
# CAFs have many genes the same as Fibro(Hughes2020)
# melanoma have many genes the same as melanocyt(Hughes2020)

## combine two sets####

se <- which(myset1$CellType %in% c("CAFs","melanoma","Macrophages",
                                   "B cells", "T cells"))
geneset <- rbind(myset1[se,],geneset)

table(geneset$CellType)

##Change cell type name
ref <- read.csv(myinf3,header = T)
refxx <- filter(ref,Dataset == "Hughes2020"|Dataset == "Tirosh2016")
refxx <- refxx[,c("Method.Name", "Paper.Name")]
sum(!unique(geneset$CellType) %in% unique(refxx$Paper.Name)) #0
sum(!unique(refxx$Paper.Name) %in% unique(geneset$CellType)) #0
geneset <- merge(geneset,refxx,by.x = "CellType", by.y = "Paper.Name")

geneset <- geneset[4:2]
setnames(geneset,old = "Method.Name", new = "CellType")
dim(geneset) # 663   3

table(geneset[1])

# CAF          Endo    Fibroblast Hair follicle  Keratinocyte    Langerhans           LEC    Melanocyte      Melanoma       Schwann      Sebocyte 
# 88            26            30            22            22            19            25            26            47            27            25 
# Skin B     Skin Mast  Skin Myeloid   Skin Plasma        Skin T          TI B TI Macrophage          TI T          VSMC 
# 21            27            18            20            29            31            92            38            30 


## save #######
CellType_Hughes2020_Tirosh2016_melanoma_TME <- geneset
# save(CellType_Hughes2020_Tirosh2016_melanoma_TME, file = myoutf1)