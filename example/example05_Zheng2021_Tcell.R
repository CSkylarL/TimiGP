#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is an example how to use TimiGP with Zheng2021 Tcell annotation
# Date: 03/07/2022
# Chenyang Skylar Li
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(TimiGP)
rm(list=ls())
#A Preprocess data #############################################################
rm(list=ls())
#1. Load SCKCM06 data ----------------------------------------------------------
data("SKCM06info")
head(SKCM06info)
data("SKCM06rna")
#2. Load cell type and marker annotation ---------------------------------------
data("CellType_Zheng2021_Tcell")
geneset <- CellType_Zheng2021_Tcell
marker <- unique(geneset$Gene)
#3. Preprocess: TimiCheckEvent & TimiPrePropress
dim(SKCM06info)
info <- TimiCheckEvent(SKCM06info)
dim(info)
dim(SKCM06rna)

rna <- TimiPrePropress(marker = marker,rna = SKCM06rna,cohort = rownames(info))
#4. Generate marker pair score: TimiGenePair  ----------------------------------
mps <- TimiGenePair(rna)
dim(mps)
#5. Perform univariate Cox regression: TimiCOX ---------------------------------

res <- TimiCOX(mps = mps,info = info,p.adj = "BH",parallel = T, core=20)
mps <- res$mps
cox_res <- res$cox_res

Tcell_COX_MP_SKCM06 <- cox_res
Tcell_MPS_SKCM06 <- mps


#save(Tcell_COX_MP_SKCM06, file = "~/Mypackage/MSofTimiGP/Fig6/Zheng2021_Tcell/Tcell_COX_MP_SKCM06.rda")
# The file is saved to the manuscript code repo: https://github.com/CSkylarL/MSofTimiGP
# You can find it here: 
# https://github.com/CSkylarL/MSofTimiGP/blob/master/Fig6/Zheng2021_Tcell/Tcell_COX_MP_SKCM06.rda


#B Gene Pair and gene network ##################################################
rm(list=ls())
# 6. Generate Directed Gene Network:TimiGeneNetwork  ----
data(Tcell_COX_MP_SKCM06)
cox_res <- Tcell_COX_MP_SKCM06

# You can use Cytoscape to visualize the network
NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Zheng2021", 
                       export = F)
#                      export =TRUE, path = "./")
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)


#B Cell interaction and network ################################################
rm(list=ls())
# 7. Generate Cell Interaction Annotation: TimiCellPair ------------------------

data("CellType_Zheng2021_Tcell")
geneset <- CellType_Zheng2021_Tcell
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

# 8. Select marker pairs A_B=1 associated with good prognosis ------------------
data(Tcell_COX_MP_SKCM06)
cox_res <- Tcell_COX_MP_SKCM06
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]

# 9. generate background: TimiBG ----
background <- TimiBG(marker.pair = row.names(cox_res))

# 10. Enrichment Analysis: TimiEnrich ----
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)

# 11. Generate Directed Cell Network:TimiCellNetwork  ----
# You can use Cytoscape to visualize the network
NET <- TimiCellNetwork(resdata = res,dataset = "Zheng2021",
                       export = F)
#                      export =TRUE, path = "./")

# 12. Visualization: Dot plot of selected cell interaction: TimiDotplot---------

p <- TimiDotplot(resdata = res,select = c(1:10))
p

# 13. Visualization: Chord Diagram of inter-cell interaction: TimiCellChord-----


# Cell Chord Diagram
TimiCellChord(resdata = res,dataset = "Zheng2021")
# Chord Diagram of marker pairs in seltect cell interaction
TimiGeneChord(resdata = res,select = 1)

# 134 Calculate favorability score: TimiFS -------------------------------------
# Visualization: TimiFSBar
# Calculate
score <- TimiFS(res)
head(score)
# Visualization
p <- TimiFSBar(score)
p
