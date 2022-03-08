# Use the immune markers which combined Charoentong2017_Galon2013_TIP
# add classic markers
rm(list = ls())
## setwd("./") # working directory is the TimiGP folder
myinf1 <- "data/CellType_Charoentong2017_Galon2013_TIP_Immune.rda"
myoutf1 <- "data/Immune_Marker_n1326.rda"
load(myinf1)

mygen <- unique(c(myset1$Gene, myset2_wo_cancer$Gene, myset3$Gene))

length(mygen) # 1302

## Add some immue gene------------------
xx <- icg <- c("ADORA2A", "CD276", "VTCN1", "BTLA", "CTLA4", "IDO1", "KIR", "LAG3", "CYBB", "PDCD1", "CD274", "HAVCR2", "SIGLEC7", "VISTA", "VSIR", "C10orf54")
mygen <- unique(c(mygen, xx))
length(mygen)			# 1307

xx <- imm.inh <- c("CTLA4", "PDCD1", "LAG3", "BTLA", "CD160", "IDO1", "IL10", "IL10RB", "TGFB1", "TGFBR1", "VTCN1", "CD244", "LGALS9", "HAVCR2", "ADORA2A", "TIGIT", "CSF1R", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KDR", "CD96", "PVRL2", "C10orf54")
mygen <- unique(c(mygen, xx))
length(mygen)			# 1310

xx <- imm.sti <- c("MICA", "MICB", "CD27", "CD274", "CD28", "CD40", "CD40LG", "CD70", "CD80", "CD86", "ICOS", "ICOSLG", "IL6", "IL6R", "PDCD1LG2", "TMEM173", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF13", "TNFSF13B", "TNFSF18", "TNFSF4", "TNFSF9", "TNFSF15", "TNFRSF25", "HHLA2", "TMIGD2", "BTNL2", "CD276", "CD48", "TNFSF14", "TNFRSF8", "PVR", "LTA",  "IL2RA", "ENTPD1", "NT5E", "CXCR4", "CXCL12", "KLRK1", "NKG2A", "RAET1E", "ULBP1")
mygen <- unique(c(mygen, xx))
length(mygen)			# 1326

xx <- imm.oth <- c("GZMA", "PRF1")
mygen <- unique(c(mygen, xx))
length(mygen)			## 1326 genes

Immune_Marker_n1326 <- mygen
# save(Immune_Marker_n1326,file = myoutf1)