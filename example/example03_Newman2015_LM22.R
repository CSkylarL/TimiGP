#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is an example how to use TimiGP with Newman2015_LM22 annotation
# Date: 03/07/2022
# Chenyang Skylar Li
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(TimiGP)
rm(list=ls())
#A Preprocess data ############################################################
rm(list=ls())
#1. Load SCKCM06 data ----------------------------------------------------------
data("SKCM06info")
head(SKCM06info)
data("SKCM06rna")
#2. Load cell type and marker annotation ---------------------------------------
data("CellType_Newman2015_LM22")
geneset <- CellType_Newman2015_LM22
marker <- unique(geneset$Gene)
#3. Preprocess: TimiCheckEvent & TimiPrePropress--------------------------------
dim(SKCM06info)
info <- TimiCheckEvent(SKCM06info)
dim(info)
dim(SKCM06rna)

rna <- TimiPrePropress(marker = marker,rna = SKCM06rna,cohort = rownames(info))
#4. Generate marker pair score: TimiGenePair  ----------------------------------
mps <- TimiGenePair(rna)
dim(mps)
#5. Perform univariate Cox regression: TimiCOX ---------------------------------

res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
mps <- res$mps
cox_res <- res$cox_res

Newman2015_COX_MP_SKCM06 <- cox_res
Newman2015_MPS_SKCM06 <- mps

# This step takes about 20-30 min, the result has been saved in data as examples
#save(Newman2015_COX_MP_SKCM06, file = "data/Newman2015_COX_MP_SKCM06.rda")



#B Gene Pair and gene network ##################################################
rm(list=ls())
# 6. Generate Directed Gene Network:TimiGeneNetwork  ---------------------------
data(Newman2015_COX_MP_SKCM06)
cox_res <- Newman2015_COX_MP_SKCM06

# You can use Cytoscape to visualize the network
NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Newman2015",
                       export = F)
#                       export =TRUE, path = "./")
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)


#B Cell interaction and network ################################################
rm(list=ls())
# 7. Generate Cell interaction Annotation: TimiCellPair ------------------------

data(CellType_Newman2015_LM22)
geneset <- CellType_Newman2015_LM22
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

# 8. Select marker pairs A_B=1 associated with good prognosis ------------------
data(Newman2015_COX_MP_SKCM06)
cox_res <- Newman2015_COX_MP_SKCM06
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]

# 9. generate background: TimiBG -----------------------------------------------
background <- TimiBG(marker.pair = row.names(cox_res))

# 10. Enrichment Analysis: TimiEnrich ------------------------------------------
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)

# 11. Generate Directed Cell Network:TimiCellNetwork  --------------------------
# You can use Cytoscape to visualize the network
NET <- TimiCellNetwork(resdata = res,dataset = "Newman2015",
                       export = F)
#                       export =TRUE, path = "./")
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)

# C Visualization###############################################################

# 12. Visualization: Dot plot of selected cell interaction: TimiDotplot---------
p <- TimiDotplot(resdata = res,select = c(1:10))
p

# 13. Visualization: Chord Diagram of inter-cell interaction: TimiCellChord-----

# Cell Chord Diagram
TimiCellChord(resdata = res,dataset = "Newman2015")
# Chord Diagram of marker pairs in seltect cell interaction
TimiGeneChord(resdata = res,select = 1)

# 14. Calculate favorability score: TimiFS -------------------------------------
# Visualization: TimiFSBar
# Calculate
score <- TimiFS(res)
head(score)
# Visualization
p <- TimiFSBar(score)
p
p1 <- TimiFSBar(score,select = c(1:5,(nrow(score)-2):nrow(score)))
p1
