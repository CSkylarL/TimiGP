#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is an example how to use TimiGP 
# with Charoentong2017_Bindea2013_Xu2018 Immune annotation
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
data(CellType_Charoentong2017_Bindea2013_Xu2018_Immune)
geneset <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune

#3. Preprocess: TimiCheckEvent & TimiPrePropress--------------------------------
dim(SKCM06info)
info <- TimiCheckEvent(SKCM06info)
dim(info)
dim(SKCM06rna)

t <- table(geneset$Dataset)
t
mps <- list()
cox_res <- list()
for (i in 1:3) {
  cat("\n",names(t)[i])
  se <- which(geneset$Dataset == names(t)[i])
  marker <- unique(geneset$Gene[se])
  rna <- TimiPrePropress(marker = marker,rna = SKCM06rna,
                         cohort = rownames(info))
  #4. Generate marker pair score: TimiGenePair  --------------------------------
  mps_tmp <- TimiGenePair(rna)
  dim(mps)
  #5. Perform univariate Cox regression: TimiCOX -------------------------------

  res<- TimiCOX(mps = mps_tmp,info = info,p.adj = "BH",parallel = T, core=20)
  mps[[ names(t)[i] ]]  <- res$mps
  cox_res[[names(t)[i]]]  <- res$cox_res
}
Immune3_COX_MP_SKCM06 <- cox_res
Immune3_MPS_SKCM06 <- mps

# save(Immune3_COX_MP_SKCM06, file = "~/Mypackage/MSofTimiGP/Fig4/Immune3_COX_MP_SKCM06.rda")
# The file is saved to the manuscript code repo: https://github.com/CSkylarL/MSofTimiGP
# You can find it here: 
# https://github.com/CSkylarL/MSofTimiGP/blob/master/Fig4/Immune3_COX_MP_SKCM06.rda


#B Gene Pair and gene network ##################################################
rm(list=ls())
# 6. Generate Directed Gene Network:TimiGeneNetwork  ---------------------------
data(Immune3_COX_MP_SKCM06)
cox_res <- Immune3_COX_MP_SKCM06
t <- names(cox_res)
# You can use Cytoscape to visualize the network
NET <- list()
for (i in 1:3){
  NET <- TimiGeneNetwork(resdata = cox_res[[t[i]]],dataset = "Immune3",
                         export = F)
  #                       export =TRUE, path = "./")
}


#B Cell interaction and network ################################################
rm(list=ls())

data("CellType_Charoentong2017_Bindea2013_Xu2018_Immune")
geneset <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune

# There are 3 genesets, find cell interaction seperately
data(Immune3_COX_MP_SKCM06)
cox_res <- Immune3_COX_MP_SKCM06
t <- names(cox_res)
res <- list()
for (i in 1:3) {
  cat("\n",t[i])
  # 7. Select marker pairs A_B=1 associated with good prognosis ----------------
  GP <- rownames(cox_res[[t[i]]])[which(cox_res[[t[i]]]$QV<0.05)]
  
  # 8. generate background: TimiBG ----
  background <- TimiBG(marker.pair = row.names(cox_res[[t[i]]]))
  
  
  # 9. Generate Cell Interaction Annotation: TimiCellPair ----------------------
  # “CellType_Charoentong2017_Bindea2013_Xu2018_Immune” 
  # contains three different immune cell-type annotations.
  # so we used "dataset" parameter to specify the dataset name.
  cell_pair <- TimiCellPair(geneset = geneset, dataset= t[i],core = 20)
  # 10. Enrichment Analysis: TimiEnrich ----
  res[[t[i]]] <- TimiEnrich(gene = GP, background = background, 
                    geneset = cell_pair, p.adj = "BH",core=20)
}


# 11. Generate Directed Cell Network:TimiCellNetwork  --------------------------
# You can use Cytoscape to visualize the network
NET <- list()
for (i in 1:3){
  NET <- TimiCellNetwork(resdata = res[[t[i]]],dataset = t[i],export = F)
  #                       export =TRUE, path = "./")
}


# 12. Visualization: Dot plot of selected cell interactions: TimiDotplot-------

p1 <- TimiDotplot(resdata = res$Charoentong2017,select = c(1:10))
p1

p2 <- TimiDotplot(resdata = res$Bindea2013,select = c(1:10))
p2

p3 <- TimiDotplot(resdata = res$Xu2018,select = c(1:10))
p3
# 13. Visualization: Chord Diagram of inter-cell interaction: TimiCellChord-----

# Cell Chord Diagram
TimiCellChord(resdata = res$Charoentong2017,dataset = "Charoentong2017")
TimiCellChord(resdata = res$Bindea2013,dataset = "Bindea2013")
TimiCellChord(resdata = res$Xu2018,dataset = "Xu2018")
# Chord Diagram of marker pairs in seltect cell interaction
TimiGeneChord(resdata = res$Charoentong2017,select = 1)
TimiGeneChord(resdata = res$Bindea2013,select = 1)
TimiGeneChord(resdata = res$Xu2018,select = 1)
# 14. Calculate favorability score: TimiFS -------------------------------------
# Visualization: TimiFSBar
# Calculate
score1 <- TimiFS(res$Charoentong2017)
score2 <- TimiFS(res$Bindea2013)
score3 <- TimiFS(res$Xu2018)
head(score)
# Visualization
p1 <- TimiFSBar(score1)
p1

p2 <- TimiFSBar(score2)
p2

p3 <- TimiFSBar(score3)
p3
