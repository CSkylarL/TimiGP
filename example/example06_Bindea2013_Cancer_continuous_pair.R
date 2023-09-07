#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is an example how to use TimiGP to capture continuous relations 
# between gene pairs with Bindea2013_cancer annotation,
# which is more sensitive compared to default logical relation.
# Date: 04/03/2023
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
data("CellType_Bindea2013_cancer")
geneset <- CellType_Bindea2013_cancer
marker <- unique(geneset$Gene)
#3. Preprocess: TimiCheckEvent & TimiPrePropress--------------------------------
dim(SKCM06info)
info <- TimiCheckEvent(SKCM06info)
dim(info)
dim(SKCM06rna)

# Note: set `GMNorm = F` to avoid default gene-wise median normalization
rna <- TimiPrePropress(marker = marker,rna = SKCM06rna,
                       cohort = rownames(info),log = T,GMNorm = F)

#4. Generate marker pair score: TimiGenePair  ----------------------------------
# Note: set `cont = T` to to capture continuous relation (expression ratio)
mps <- TimiGenePair(rna, cont = T) 
dim(mps)
#5. Perform univariate Cox regression: TimiCOX ---------------------------------

res <- TimiCOX(mps = mps,info = info,p.adj = "BH",parallel = T, core=20)
mps <- res$mps
cox_res <- res$cox_res

Bindea2013c_COX_MP_SKCM06_CONT <- cox_res
Bindea2013c_MPS_SKCM06_CONT <- mps

# This step takes about 5-10 min, please be patient


#B Gene Pair and gene network ##################################################
# 6. Generate Directed Gene Network:TimiGeneNetwork  ---------------------------

# You can use Cytoscape to visualize the network
NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Bindea2013_Cancer",
                       export = F)
#                       export =TRUE, path = "./")
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)


#B Cell interaction and network ################################################
# 7. Generate Cell interaction Annotation: TimiCellPair ------------------------

data(CellType_Bindea2013_cancer)
geneset <- CellType_Bindea2013_cancer
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

# 8. Select marker pairs A_B=1 associated with good prognosis ------------------
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]

# 9. generate background: TimiBG -----------------------------------------------
background <- TimiBG(marker.pair = row.names(cox_res))

# 10. Enrichment Analysis: TimiEnrich ------------------------------------------
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)

# 11. Generate Directed Cell Network:TimiCellNetwork  --------------------------
# You can use Cytoscape to visualize the network

NET <- TimiCellNetwork(resdata = res,dataset = "Bindea2013_Cancer",
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
TimiCellChord(resdata = res,dataset = "Bindea2013_Cancer")
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

