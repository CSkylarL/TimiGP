# TimiGP
 A R package to infer molecular and cellular interactions in tumor immune microenvironment through gene pairs
## Installation
 1. Install [the `devtools` package](https://github.com/r-lib/devtools) from CRAN.
 2. Run `devtools::install_github("CSkylarL/TimiGP")`.
 
 Follow below code:
```R
# Install devtools from CRAN.
install.packages("devtools")

# Install compassR from YosefLab.
devtools::install_github("CSkylarL/TimiGP")
```
## Usage
Here is a summary of the functions in the package:
| Step                                                               | Function                                                        |
| ------------------------------------------------------------------ | --------------------------------------------------------------- |
| 1.1   Preprocess of Clinical Info                                  | `TimiCheckEvent`(info = NULL) |  
| 1.2   Preprocess of Transcriptome Profile                          | `TimiPrePropress`(marker, rna = NULL, cohort, log = TRUE)|
| 2.1   Pairwise Comparison                                          | `TimiGenePair`(rna = NULL)|
| 2.2   Directed IRGP Selection                                      | `TimiCOX`(mps = NULL,info = NULL, p.adj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))|
| 2.3   Directed Gene Network                                        | `TimiGeneNetwork`(resdata = NULL, select = NULL, dataset = NULL, geneset = NULL, export = TRUE, path = NULL)|
| 3.1   Cell Pair Annotation                                         | `TimiCellPair`(geneset = NULL, dataset = NULL, core = 1)|
| 3.2   Prepare Enrichment Background                                | `TimiBG`(marker.pair = NULL)|
| 3.3.1 Enrichment Analysis                                          | `TimiEnrich`(gene = NULL, background = NULL, geneset = NULL, p.adj = NULL, core = 1, pair = TRUE)|
| 3.3.2 Visualization: Dotplot of Enrichment                         | `TimiDotplot`(resdata = NULL, select = 1:5)|
| 3.4.1 Visualization: Chord Diagram of Cell Interaction             | `TimiCellChord`(resdata = NULL,select = NULL, dataset = NULL,group = NULL, color = NULL)|
| 3.4.2 Visualization: Chord Diagram of Selected Gene Interaction    | `TimiGeneChord`(resdata = NULL, select = 1, color = NULL)|
| 3.5.1 Cell Interaction Network                                     | `TimiCellNetwork`(resdata = NULL, select = NULL, dataset = NULL, group = NULL, geneset = NULL, export = TRUE, path = NULL)|
| 3.5.2 Favorability Score                                           | `TimiFS`(resdata = NULL, cutoff = 0.05)|
| 3.5.3 Visualization: Bar plot of Favorability Score                | `TimiFSBar`(score = NULL, select = NULL)|

## Available Data
Here is a summary of available data in the package:
| Available Data                                             | Description                                              |
| ---------------------------------------------------------- | -------------------------------------------------------- | 
| 1.1   data(SKCM06info)                                     | TCGA SKCM06(metastatic melnoma) clinical info                                                                             |
| 1.2   data(SKCM06rna)                                      | TCGA SKCM06(metastatic melnoma) transcriptome profile                                                                     |
| 2.1   data(Immune_Marker_n1326)                            | #1326 immune cell markers                                                                                                 |
| 2.2   data(`CellType_Charoentong2017_Galon2013_TIP_Immune`)| Immune cell types and markers                                                                                             |
| 2.3   data(Immune3_COX_MP_SKCM06)                          | `TimiCOX` results with above annotation that reveals the association between each marker pairs and favorable prognosis.   |
| 3.1   data(`CellType_Alex2020_Levi2019_TME`)               | Melanoma tumor microenvironment cell types and markers                                                                    |
| 3.2   data(TME_COX_MP_SKCM06)                              | `TimiCOX` results with above annotation that reveals the association between each marker pairs and favorable prognosis.   |
| 4.1   data(`CellType_Zhang2021S2_Tcell`)                   | T cell subtypes and markers                                                                                               |
| 4.2   data(Tcell_COX_MP_SKCM06)                            | `TimiCOX` results with above annotation that reveals the association between each marker pairs and favorable prognosis.   |
| 5.1   data(`CellType_Galon2013_cancer`)                    | Markers of immune cell types and 1 cancer cell type                                                                       |
| 5.2   data(Galon2013c_COX_MP_SKCM06)                       | `TimiCOX` results with above annotation that reveals the association between each marker pairs and favorable prognosis.   |
| 5.3   data(Galon2013c_enrich)                              | An example of `TimiEnrich` results on cell interactions                                                                   |

## Example
Here is an example that how to use the package to infer gene and cell interaction based on relative abundance. 
Other examples can be found in the `example` folder. And the process to generate the cell markers can be found in `inst/extdata`.

Library the package
```R
library(TimiGP)
rm(list=ls())
```

### 1.   Preprocess of clinical infomation and transctiptome profile
`TimiCheckEvent` will filter out the low quality clincal info(event and time-to-event).

`TimiPrePropress` will extracted the transcription profile of selected markers and cohorts, then perform log transformation and gene wise median normalization.
```R
rm(list=ls())
#1. Load SCKCM06 data ----
data("SKCM06info")
head(SKCM06info)
data("SKCM06rna")
#2. Load cell type and marker annotation ----
data("CellType_Galon2013_cancer")
geneset <- CellType_Galon2013_cancer
marker <- unique(geneset$Gene)
#3. Preprocess: TimiCheckEvent & TimiPrePropress ----
info <- TimiCheckEvent(SKCM06info)
rna <- TimiPrePropress(marker = marker,rna = SKCM06rna,cohort = rownames(info))
```
### 2.1   Pairwise Comparison
`TimiGenePair` will Capture logical relation of any two marker pairs, and generate a matrix of Marker Pair Score(MPS):

  - 1 or TRUE = the expression of gene A > that of gene B, 

  - 0 or FALSE = the expression of gene A < that of gene B.
```R
#4. Generate marker pair score: TimiGenePair  ----
mps <- TimiGenePair(rna)
#5. Perform univariate Cox regression to find the association between marker pair and survival: TimiCOX ----
```
res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
mps <- res[[1]]
cox_res <- res[[2]]

Galon2013c_COX_MP_SKCM06 <- cox_res
Galon2013c_MPS_SKCM06 <- mps

# This step takes about 20-30 min, the result has been saved in data as examples
#save(Galon2013c_COX_MP_SKCM06, file = "data/Galon2013c_COX_MP_SKCM06.rda")


xx <- cox_res[order(cox_res[,2]),]
sum(xx[,3]<0.05)  # 2644
sum(xx[,3]<0.05)/nrow(xx) # 2%
remove(xx)



#B Gene Pair and gene network ###################################################
rm(list=ls())
# 6. Generate Directed Gene Network:TimiGeneNetwork  ----
data(Galon2013c_COX_MP_SKCM06)
cox_res <- Galon2013c_COX_MP_SKCM06
# I saved the output files in notebook. Please choose yours
# You can use Cytoscape to visualize the network
NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Galon2013_Cancer",export =TRUE, path = "./")
head(NET$network,n = 3)
head(NET$node,n = 3)
head(NET$edge,n = 3)


#B Cell interaction and network ###########################################
rm(list=ls())
# 7. Generate Cell Pair Annotation: TimiCellPair ----

data(CellType_Galon2013_cancer)
geneset <- CellType_Galon2013_cancer
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

# 8. Select marker pairs A_B=1 associated with good prognosis ----
data(Galon2013c_COX_MP_SKCM06)
cox_res <- Galon2013c_COX_MP_SKCM06
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]

# 9. generate background: TimiBG ----
background <- TimiBG(marker.pair = row.names(cox_res))

# 10. Enrichment Analysis: TimiEnrich ----
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)
# This has been saved to data as an example
#Galon2013c_enrich <- res
#save(Galon2013c_enrich,file = "data/Galon2013c_enrich.rda")

# 11. Generate Directed Cell Network:TimiCellNetwork  ----
# You can use Cytoscape to visualize the network
rm(list=ls())
data("Galon2013c_enrich")
res <- Galon2013c_enrich
NET <- TimiCellNetwork(resdata = res,dataset = "Galon2013_Cancer",export =TRUE, path = "./")

# C Visualization##########################################################

# 12. Visualization: Dot plot of selected cell pair enrichment: TimiDotplot-----
rm(list=ls())
data("Galon2013c_enrich")
res <- Galon2013c_enrich
p <- TimiDotplot(resdata = res,select = c(1:10))
p

# 13. Visualization: Chord Diagram of significant cell pair enrichment: TimiCellChord----
rm(list=ls())
data("Galon2013c_enrich")
res <- Galon2013c_enrich
# Cell Chord Diagram
TimiCellChord(resdata = res,dataset = "Galon2013_Cancer")
# Chord Diagram of marker pairs in seltect cell pair
TimiGeneChord(resdata = res,select = 1)

# 13. Calculate favorability score: TimiFS ----
# Visualization: Timi 
rm(list=ls())
data("Galon2013c_enrich")
res <- Galon2013c_enrich
# Calculate
score <- TimiFS(res)
head(score)
# Visualization
p <- TimiFSBar(score)
p

