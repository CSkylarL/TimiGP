# Orgnize immune related genes
library(dplyr)
rm(list = ls())
myinf1 <- "inst/extdata/Charoentong2017_tableS6_Immune.txt"
myinf2 <- "inst/extdata/Bindea2013_tableS1_immune.txt"
myinf3 <- "inst/extdata/Xu2018_Signature.txt"
myoutf1 <- "data/Immune_Marker_n1293.rda"


## Charoentong2017 ----
myset1 <- read.table(myinf1, 
                     header = T, sep = "\t", stringsAsFactors = F)
head(myset1)
myset1 <- myset1[, -3]
colnames(myset1) <- c("Gene", "CellType")
myset1$Dataset <- "Charoentong2017"
myset1 <- myset1[, 3:1]

## Bindea2013 ----
myset2 <- read.table(myinf2, 
                     header = F, sep = "\t", stringsAsFactors = F)
head(myset2)
myset2 <- myset2[, -3]
colnames(myset2) <- c("CellType", "Gene")
myset2$Dataset <- "Bindea2013"
myset2 <- myset2[, c(3, 1,2)]

se <- which(myset2$CellType == "SW480 cancer cells" | 
              myset2$CellType == "Lymph vessels" |
              myset2$CellType == "Normal mucosa" |
              myset2$CellType == "Blood vessels")
length(se)
myset2 <- myset2[-se, ]
## Xu2018 ----

myset3 <- read.table(myinf3, 
                     header = T, sep = "\t", stringsAsFactors = F)
head(myset3)
myset3 <- myset3[, c(1,4)]
colnames(myset3) <- c("Gene", "CellType")
myset3$Dataset <- "Xu2018"
myset3 <- myset3[, 3:1]



mygen <- unique(c(myset1$Gene,myset2$Gene,myset3$Gene))

length(mygen) # 1269

## Add some immue gene------------------
xx <- icg <- c("ADORA2A", "CD276", "VTCN1", "BTLA", "CTLA4", "IDO1", "KIR", "LAG3", "CYBB", "PDCD1", "CD274", "HAVCR2", "SIGLEC7", "VISTA", "VSIR", "C10orf54")
mygen <- unique(c(mygen, xx))
length(mygen)			# 1274

xx <- imm.inh <- c("CTLA4", "PDCD1", "LAG3", "BTLA", "CD160", "IDO1", "IL10", "IL10RB", "TGFB1", "TGFBR1", "VTCN1", "CD244", "LGALS9", "HAVCR2", "ADORA2A", "TIGIT", "CSF1R", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KDR", "CD96", "PVRL2", "C10orf54")
mygen <- unique(c(mygen, xx))
length(mygen)			# 1277

xx <- imm.sti <- c("MICA", "MICB", "CD27", "CD274", "CD28", "CD40", "CD40LG", "CD70", "CD80", "CD86", "ICOS", "ICOSLG", "IL6", "IL6R", "PDCD1LG2", "TMEM173", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF13", "TNFSF13B", "TNFSF18", "TNFSF4", "TNFSF9", "TNFSF15", "TNFRSF25", "HHLA2", "TMIGD2", "BTNL2", "CD276", "CD48", "TNFSF14", "TNFRSF8", "PVR", "LTA",  "IL2RA", "ENTPD1", "NT5E", "CXCR4", "CXCL12", "KLRK1", "NKG2A", "RAET1E", "ULBP1")
mygen <- unique(c(mygen, xx))
length(mygen)			# 1293

xx <- imm.oth <- c("GZMA", "PRF1")
mygen <- unique(c(mygen, xx))
length(mygen)			## 1293 genes

Immune_Marker_n1293 <- mygen
# save(Immune_Marker_n1293,file = myoutf1)

