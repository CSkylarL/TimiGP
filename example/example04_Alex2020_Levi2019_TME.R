#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is an example how to use TimiGP with Alex2020_Levi2019 TME annotation
# Date: 03/07/2022
# Chenyang Li
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(TimiGP)
rm(list=ls())
#A Preprocess data #####################################################
rm(list=ls())
#1. Load SCKCM06 data ----
data("SKCM06info")
head(SKCM06info)
data("SKCM06rna")
#2. Load cell type and marker annotation ----
data("CellType_Alex2020_Levi2019_TME")
geneset <- CellType_Alex2020_Levi2019_TME
marker <- unique(geneset$Gene)
#3. Preprocess: TimiCheckEvent & TimiPrePropress
dim(SKCM06info)
info <- TimiCheckEvent(SKCM06info)
dim(info)
dim(SKCM06rna)

rna <- TimiPrePropress(marker = marker,rna = SKCM06rna,cohort = rownames(info))
#4. Generate marker pair score: TimiGenePair  ----
mps <- TimiGenePair(rna)
dim(mps)
#5. Perform univariate Cox regression to find the association between marker pair and survival: TimiCOX ----

res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
mps <- res[[1]]
cox_res <- res[[2]]

TME_COX_MP_SKCM06 <- cox_res
TME_MPS_SKCM06 <- mps

# This step takes about 20-30 min, the result has been saved in data as examples
#save(TME_COX_MP_SKCM06, file = "data/TME_COX_MP_SKCM06.rda")

xx <- cox_res[order(cox_res[,2]),]
sum(xx[,3]<0.05)  # 12742
sum(xx[,3]<0.05)/nrow(xx) # 8%
remove(xx)



#B Gene Pair and gene network ###################################################
rm(list=ls())
# 6. Generate Directed Gene Network:TimiGeneNetwork  ----
data(TME_COX_MP_SKCM06)
cox_res <- TME_COX_MP_SKCM06
# I saved the output files in notebook. Please choose yours
# You can use Cytoscape to visualize the network
NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Galon2013_Cancer",path = "./")
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)


#B Cell interaction and network ###########################################
rm(list=ls())
# 7. Generate Cell Pair Annotation: TimiCellPair ----

data("CellType_Alex2020_Levi2019_TME")
geneset <- CellType_Alex2020_Levi2019_TME
# Rename the dataset to make it one set for cell pair
geneset$Dataset <- "Alex2020_Levi2019"
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

# 8. Select marker pairs A_B=1 associated with good prognosis ----
data(TME_COX_MP_SKCM06)
cox_res <- TME_COX_MP_SKCM06
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]

# 9. generate background: TimiBG ----
background <- TimiBG(marker.pair = row.names(cox_res))

# 10. Enrichment Analysis: TimiEnrich ----
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)

# 11. Generate Directed Cell Network:TimiCellNetwork  ----
# You can use Cytoscape to visualize the network
NET <- TimiCellNetwork(resdata = res,dataset = "Alex2020_Levi2019",path = "./")

# 12. Visualization: Dot plot of selected cell pair enrichment: TimiDotplot-----

p <- TimiDotplot(resdata = res,select = c(1:10))
p

# 13. Visualization: Chord Diagram of significant cell pair enrichment: TimiCellChord----

# Cell Chord Diagram
TimiCellChord(resdata = res,dataset = "Alex2020_Levi2019")
# Chord Diagram of marker pairs in seltect cell pair
TimiGeneChord(resdata = res,select = 1)

# 13. Calculate favorability score: TimiFS ----
# Visualization: Timi 
# Calculate
score <- TimiFS(res)
head(score)
# Visualization
p <- TimiFSBar(score)
p
