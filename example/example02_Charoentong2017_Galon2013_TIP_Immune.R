#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is an example how to use TimiGP with Charoentong2017_Galon2013_TIP Immune annotation
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
#2. Load immune cell marker  ----
data("Immune_Marker_n1326")
marker <-Immune_Marker_n1326
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

Immune3_COX_MP_SKCM06 <- cox_res
Immune3_MPS_SKCM06 <- mps

# This step takes about 20-30 min, the result has been saved in data as examples
#save(Immune3_COX_MP_SKCM06, file = "data/Immune3_COX_MP_SKCM06.rda")

xx <- cox_res[order(cox_res[,2]),]
sum(xx[,3]<0.05)  # 37665
sum(xx[,3]<0.05)/nrow(xx) # 5%
remove(xx)



#B Gene Pair and gene network ###################################################
rm(list=ls())
# 6. Generate Directed Gene Network:TimiGeneNetwork  ----
data(Immune3_COX_MP_SKCM06)
cox_res <- Immune3_COX_MP_SKCM06

# You can use Cytoscape to visualize the network
NET1 <- TimiGeneNetwork(resdata = cox_res,dataset = "Immune3",export =FALSE, path = "./")


#B Cell interaction and network ###########################################
rm(list=ls())

data("CellType_Charoentong2017_Galon2013_TIP_Immune")
geneset <- CellType_Charoentong2017_Galon2013_TIP_Immune

# 7. Select marker pairs A_B=1 associated with good prognosis ----
data(Immune3_COX_MP_SKCM06)
cox_res <- Immune3_COX_MP_SKCM06
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]


# 8. generate background: TimiBG ----
background <- TimiBG(marker.pair = row.names(cox_res))

# There are 3 genesets, find cell interaction seperately
t <- table(geneset$Dataset)
t

# Charoentong2017       Galon2013             TIP 
#            782             542              85 
res <- list()
for (i in 1:3) {
  # 9. Generate Cell Pair Annotation: TimiCellPair ----
  cat("\n",names(t)[i])
  se <- which(geneset$Dataset == names(t)[i])
  cell_pair <- TimiCellPair(geneset = geneset[se,],core = 20)
  
  # 10. Enrichment Analysis: TimiEnrich ----
  res[[names(t)[i]]] <- TimiEnrich(gene = GP, background = background, 
                    geneset = cell_pair, p.adj = "BH",core=20)
}


# 11. Generate Directed Cell Network:TimiCellNetwork  ----
# You can use Cytoscape to visualize the network
NET1 <- TimiCellNetwork(resdata = res$Charoentong2017,dataset = "Charoentong2017",export =FALSE, path = "./")
NET2 <- TimiCellNetwork(resdata = res$Galon2013,dataset = "Galon2013",export =FALSE, path = "./")
NET3 <- TimiCellNetwork(resdata = res$TIP,dataset = "TIP",export =FALSE, path = "./")

# 12. Visualization: Dot plot of selected cell pair enrichment: TimiDotplot-----

p1 <- TimiDotplot(resdata = res$Charoentong2017,select = c(1:10))
p1

p2 <- TimiDotplot(resdata = res$Galon2013,select = c(1:10))
p2

p3 <- TimiDotplot(resdata = res$TIP,select = c(1:10))
p3
# 13. Visualization: Chord Diagram of significant cell pair enrichment: TimiCellChord----

# Cell Chord Diagram
TimiCellChord(resdata = res$Charoentong2017,dataset = "Charoentong2017")
TimiCellChord(resdata = res$Galon2013,dataset = "Galon2013")
TimiCellChord(resdata = res$TIP,dataset = "TIP")
# Chord Diagram of marker pairs in seltect cell pair
TimiGeneChord(resdata = res$Charoentong2017,select = 1)
TimiGeneChord(resdata = res$Galon2013,select = 1)
TimiGeneChord(resdata = res$TIP,select = 1)
# 13. Calculate favorability score: TimiFS ----
# Visualization: Timi 
# Calculate
score1 <- TimiFS(res$Charoentong2017)
score2 <- TimiFS(res$Galon2013)
score3 <- TimiFS(res$TIP)
head(score)
# Visualization
p1 <- TimiFSBar(score1)
p1

p2 <- TimiFSBar(score2)
p2

p3 <- TimiFSBar(score3)
p3
