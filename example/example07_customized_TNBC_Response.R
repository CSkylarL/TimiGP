#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is an example how to use TimiGP-Response
# Date: 06/05/2023
# Author: Chenyang Skylar Li
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Clear all variables
rm(list=ls())
# load the TimiGP package
library(TimiGP)
#A Preprocess data ############################################################
#1. Load example data ----------------------------------------------------------
# This is example data was generated from GSE194040, 
# including 26 Triple-negative breast cancer (TNBC) patients
# treated with the combination of anti-PD1 immunotherapy and chemotherpay
# Please refer to ./inst/extdata/build_TNBCaPD1_RNA_info.R for more details
data(TNBCaPD1info)
head(TNBCaPD1info)
# 1=Responder; 0=Non-responder
#           Response
# X111038        0
# X115638        0
# X116603        0
# X153327        1
# X178649        1
# X196024        0

data(TNBCaPD1rna)
#2. Load cell type and marker annotation ---------------------------------------
# The cell type annotation is customized for TNBC (anti-PD1/PD-L1) 
# from single cell RNA-seq.
# For more details, please refer:
# https://github.com/CSkylarL/MSofTimiGP-Response/blob/master/Fig1/code/Fig1.script2_generate_TimiGP_markers.R

data("CellType_TNBC_aPDL1")
geneset <- CellType_TNBC_aPDL1
marker <- unique(geneset$Gene)

#3. Preprocess:  TimiPrePropress------------------------------------------------
info <- TNBCaPD1info
table(info)

# Response
# 0  1 
# 13 13 

# Automatically check if the expression has been log tranformed or not
rna <- TNBCaPD1rna
if (max(rna,na.rm = T) < 100){
  log = FALSE
} else {
  log = TRUE
}
cat("\n", max(rna,na.rm = T)," log = ", log,"\n")
# 18.565  log =  FALSE 

rna <- TimiPrePropress(marker = marker, cohort = rownames(info),
                       log = log, GMNorm = T, rna = rna)
#4. Generate marker pair score: TimiGenePair  ----------------------------------
mps <- TimiGenePair(rna)
dim(mps)
#5. Perform Fisher Test: TimiFisher --------------------------------------------
# Note: This step is unique for the TimiGP-Response
# If you have >1000 markers, 
# you can set `parallel = T, core=10` to enable parallel computiong

res <- TimiFisher(mps = mps,info = info,p.adj = "BH")
  
mps <- res$mps
fisher_res <- res$fisher_res

#B Gene Pair and gene network ##################################################
# 6. Generate Directed Gene Network:TimiGeneNetwork  ---------------------------

# You can use Cytoscape to visualize the network
NET <- TimiGeneNetwork(resdata = fisher_res,
                       dataset = "Other",
                       geneset = geneset,
                       condition = "PV",
                       cutoff = 0.05,
                       export = F)
# If you want to export the files,
# Please set `export =TRUE, path = "./"`
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)


#B Cell interaction and network ################################################
# 7. Generate Cell interaction Annotation: TimiCellPair ------------------------
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

# 8. Select marker pairs A_B=1 associated with responders ------------------

GP <- fisher_res %>% filter(PV < 0.05)  %>% rownames()

# 9. generate background: TimiBG -----------------------------------------------
background <- TimiBG(marker.pair = row.names(fisher_res))

# 10. Enrichment Analysis: TimiEnrich ------------------------------------------
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)

# Optional: Permutation to generate FDR: TimiPermFFDR --------------------------
res <- TimiPermFDR(resdata = res, geneset = geneset, gene = GP,
                   background = background, niter = 100, core = 20)

# 11. Generate Directed Cell Network:TimiCellNetwork  --------------------------
# define cell-cell interaction -------------------------------------------
CI_condition <- "Adjust.P.Value"
CI_cutoff <- 0.05
cat("CI_cutoff: ", CI_condition, " < ",
    CI_cutoff, " : #",sum(res[,CI_condition] < CI_cutoff))
# CI_cutoff:  Adjust.P.Value  <  0.05  : # 49

# settings for visualization  --------------------------------------------------
# Cell names
cell <- c("CD20+B" , "CD79a+Plasma",
          "CD4+TCF1+T" ,"CD4+PD1+T" ,"Treg" , 
          "CD8+TCF1+T" ,"CD8+GZMB+T" ,"CD8+PD1+Tex" ,
          "CD56+NK" ,"cDC" , "pDC" , 
          "M1" , "M2" , 
          "Mast" ) 

# Cell group
group <- c(rep("B",2),
           rep("CD4T", 3),
           rep("CD8T", 3),
           "NK",
           rep("DC", 2),
           rep("Mac", 2),
           "Mast")
# Cell color
color <- c("#66c2a4","#238b45",
           "#bcbddc","#807dba","#54278f",
           "#fa9fb5","#dd3497","#7a0177",
           "#2b8cbe",
           "#225ea8","#081d58",
           "#feb24c","#ec7014",
           "#e31a1c")
names(group) <- cell
names(color) <- cell

# You can use Cytoscape to visualize the network

NET <- TimiCellNetwork(resdata = res,
                       condition = CI_condition, 
                       cutoff = CI_cutoff,
                       geneset = geneset,
                       dataset = "Other",
                       group = group,
                       export =F)
# If you want to export the files,
# Please set `export =TRUE, path = "./"`
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)

# C Visualization###############################################################

# 12. Visualization: Dot plot of selected cell interaction: TimiDotplot---------
nn <- min(sum(res[,CI_condition] < CI_cutoff),10)
p <- TimiDotplot(resdata = res,
                 condition = CI_condition,
                 cutoff = CI_cutoff,
                 select = c(1:nn))
p

# 13. Visualization: Chord Diagram of inter-cell interaction: TimiCellChord-----
TimiCellChord(resdata = res,
              dataset = "Other",group = group,color = color,
              condition = CI_condition, cutoff = CI_cutoff)
# Chord Diagram of marker pairs in seltect cell interaction
TimiGeneChord(resdata = res,select = 1)

# 14. Calculate favorability score: TimiFS -------------------------------------
# Visualization: TimiFSBar
score <- TimiFS(res,condition = CI_condition,cutoff = CI_cutoff)
head(score)
# Visualization
p <- TimiFSBar(score)
p

# This is an example code modified from
# https://github.com/CSkylarL/MSofTimiGP-Response/blob/master/Fig1/code/Fig1.script3_customized_TimiGP_response_for_TNBC.R
# The detailed output can be found from 
# https://github.com/CSkylarL/MSofTimiGP-Response/tree/master/Fig1/result/TimiGP_TNBC_scRNA_marker_for_spatial/Wolf_GSE194040_ISPY2_ptx_pem_TNBC