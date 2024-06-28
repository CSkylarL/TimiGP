# TimiGP

[![DOI](https://zenodo.org/badge/467312934.svg)](https://zenodo.org/badge/latestdoi/467312934)

`TimiGP (Tumor Immune Microenvironment Illustration based on Gene Pairing)` is a computational framework to infer cell-cell interactions and clinical values in tumor immune microenvironment through gene pairs. This repository contains the R package of TimiGP, which is intended for research use only.

For more details, please jump to the [Publication & Resources](#publication--resources) section.

## Table of Contents
- [Rationale](#rationale)
- [Framework](#framework)
- [Application](#application)
- [Citation](#citation)
- [Publication & Resources](#publication--resources)
- [TimiGP system To-do-list](#timigp-system-to-do-list)
- [Contributing](#contributing)
- [Installation](#installation)
- [Function](#function)
- [Available Data](#available-data)
- [Example](#example)

## Rationale

The tumor immune microenvironment(TIME) is a balance between anti-tumor and pro-tumor immune cells. If the function of anti-tumor cell types is more vital than the pro-tumor cells (e.g., higher abundance, higher marker expression), the TIME is associated with favorable patients’ clinical outcome; otherwise, it is associated with unfavorable outcome.

![rationale](/assets/images/rationale.png)

## Framework
Here is the **Overview of TimiGP framework**.

![Overview of TimiGP framework](/assets/images/framework.png)

`Inputs`
 1. Cell-type marker (single cell RNA-seq, prior knowledge, etc.);
 2. Bulk Transcriptomic Profile (RNA-seq or microarray);
 3. Clinical informatoin:
    - `Prognosis Module`: `Survival statistics`, including events (e.g., death, recurrence) and time-to-event of the same cohort. 
    - `Response Module`: `Therapy response`, including response (e.g., 1=respond, 0=non-respond) of the same cohort.
    - `Single Sample Module`: `none` (It has no been released yet, please stay tuned).

`TimiGP Steps`
 1. Define and select marker gene pair matrix;
 2. Choose gene pairs associated with favorable prognosis (Prognosis Module) or response to immunotherapy (Response Module);
 3. Construct the directed gene-gene network; 
 4. Identify cell-cell interactions by enrichment analysis;
 5. Build and analyze the cell-cell interaction network.
 
`Outputs`
 1. Cell-cell interaction network,
 2. Clinical value of each cell associated with the input outcome (favorability score).
 
## Application
TimiGP is designed to **infer the cell-cell interaction network and clinical association of immune cells**. Based on the resulting immunological insights, The method will **facilitate the development of prognostic/predictive models**. Taking advantage of different biomarker references derived from **bulk and single-cell RNA-seq**, TimiGP can be applied to investigate the **entire tumor microenvironment** or **cell subpopulations**, perform **pan-cancer analysis** or study **other diseases**. 

 ![application](/assets/images/application.png)
 
## Citation
This package is intended for research use only. 

If you use TimiGP-Prognosis in your publication, please cite the paper:

Li, C. et al. TimiGP: Inferring cell-cell interactions and prognostic associations in the tumor immune microenvironment through gene pairs. Cell Rep Med 4, 101121, doi:10.1016/j.xcrm.2023.101121 (2023).

If you use TimiGP-Response in your publication, please cite the paper:

Li, C. et al. TimiGP-Response: the pan-cancer immune landscape associated with response to immunotherapy. bioRxiv, 2024.2006.2021.600089, doi:10.1101/2024.06.21.600089 (2024).

For the detailed protocol, please refer to:

Li, C., Zhang, J. & Cheng, C. TimiGP: An R package to depict the tumor microenvironment from bulk transcriptomics. STAR protocols 4, 102742 (2023).

## Publication & Resources
`TimiGP v1.1.0 Prognosis Module:` Base version
- Cell Reports Medicine: [TimiGP: Inferring cell-cell interactions and prognostic associations in the tumor immune microenvironment through gene pairs](https://doi.org/10.1016/j.xcrm.2023.101121)
- Code for the above paper: [MSofTimiGP](https://github.com/CSkylarL/MSofTimiGP)
- bioRxiv preprint: [TimiGP: inferring cell-cell functional interactions and clinical values in the tumor immune microenvironment through gene pairs.](https://www.biorxiv.org/content/10.1101/2022.11.17.515465v1.full)

`TimiGP v1.2.0 Prognosis Module:` Optimized computational efficiency
- STAR Protocols: [TimiGP: An R package to depict the tumor microenvironment from bulk transcriptomic](https://star-protocols.cell.com/protocols/3156)

`TimiGP v1.3.0 Response Module`:
- bioRxiv preprint: [TimiGP-Response: the pan-cancer immune landscape associated with response to immunotherapy](https://www.biorxiv.org/content/10.1101/2024.06.21.600089v1.full)
- Code & Resources: [MSofTimiGP-Response](https://github.com/CSkylarL/MSofTimiGP-Response)
- ASHG 2023 abstract (Page 135): [Utilizing TimiGP for in-depth analysis of the tumor immune microenvironment and its association with clinical outcomes of various cancers after immunotherapy](https://www.ashg.org/wp-content/uploads/2023/10/ASHG2023-PlatformAbstracts.pdf)
- AACR 2023 abstract: [TimiGP: Dissect the tumor immune microenvironment and its association with survival and immunotherapy response](https://aacrjournals.org/cancerres/article/83/7_Supplement/2080/723394/Abstract-2080-TimiGP-Dissect-the-tumor-immune)

`TimiGP v1.4.0 Single Sample Module`:
- AACR 2024 abstract: [Equip the TimiGP system with a single-sample module, facilitating the cell-cell interaction inference of the tumor microenvironment in any sample](https://doi.org/10.1158/1538-7445.AM2024-4946)

## TimiGP system To-do-list

- [x] TimiGP concept: [prognosis module](https://doi.org/10.1016/j.xcrm.2023.101121) and R package (v1.1.0)
- [x] Code optimization (v1.2.0) and publish the [protocol for the prognosis module](https://star-protocols.cell.com/protocols/3156)
- [x] TimiGP analysis: response module (wrap-up)
- [x] TimiGP tool: **single sample module** (validated)
- [X] Publish preprint and package for **response module** 
- [ ] Publish peer-reviewed paper of **response module** 
- [ ] TimiGP upgrade: **single cell module** (conceptualized)
- [ ] TimiGP WebApp

## Contributing

We greatly welcome contributions to TimiGP. Please submit a pull request if you have any ideas (e.g., new modules and functions) or bug fixes. We also welcome any issues you encounter while using TimiGP.

## Installation
 1. Install [the `devtools` package](https://github.com/r-lib/devtools) from CRAN.
 2. Run `devtools::install_github("CSkylarL/TimiGP")`.
 
 Follow the below code:
```R
# Install devtools from CRAN.
install.packages("devtools")

# Install TimiGP from GitHub.
devtools::install_github("CSkylarL/TimiGP")

# Load the package
library(TimiGP)
```
## Function
Here is a summary of the functions in the package:
| Step                                                               | Function                                                        |
| ------------------------------------------------------------------ | --------------------------------------------------------------- |
| 1.1   Preprocess of Clinical Info                                  | `TimiCheckEvent`(info = NULL) |  
| 1.2   Preprocess of Transcriptome Profile                          | `TimiPrePropress`(marker, rna = NULL, cohort, log = TRUE, GMNorm=TRUE)|
| 2.1   Pairwise Comparison                                          | `TimiGenePair`(rna = NULL, cont = FALSE)|
| 2.2.1   Prognostic IMGP Selection with Cox Regression (Prognosis Module)         | `TimiCOX`(mps = NULL,info = NULL, p.adj = "BH", parallel = FALSE, core = 1)|
| 2.2.2   Responsive IMGP Selection with Fisher Test (Response Module)             | `TimiFisher`(mps = NULL, info = NULL, p.adj = "BH", parallel = FALSE, core = 1)|
| 2.3   Directed Gene Network                                        | `TimiGeneNetwork`(resdata = NULL, select = NULL, dataset = NULL, geneset = NULL, condition = "QV", cutoff = 0.05, export = TRUE, path = NULL)|
| 3.1   Cell Interaction Annotation                                  | `TimiCellPair`(geneset = NULL, dataset = NULL, core = 1)|
| 3.2   Prepare Enrichment Background                                | `TimiBG`(marker.pair = NULL)|
| 3.3.1 Enrichment Analysis                                          | `TimiEnrich`(gene = NULL, background = NULL, geneset = NULL, p.adj = NULL, core = 1, pair = TRUE)|
| 3.3.2 Permutation Test                                             | `TimiPermFDR`(resdata = NULL, geneset = NULL, gene = NULL, background = NULL, niter = 100, core = 1)|
| 3.3.3 Visualization: Dotplot of Cell Interaction                   | `TimiDotplot`(resdata = NULL, select = 1:5)|
| 3.4.1 Visualization: Chord Diagram of cell-cell Interaction       | `TimiCellChord`(resdata = NULL,select = NULL, dataset = NULL,group = NULL, color = NULL, condition = "Adjust.P.Value", cutoff = 0.05)|
| 3.4.2 Visualization: Chord Diagram of Selected Gene Interaction    | `TimiGeneChord`(resdata = NULL, select = 1, color = NULL)|
| 3.5.1 Functional Interaction Network                               | `TimiCellNetwork`(resdata = NULL, select = NULL, dataset = NULL, group = NULL, geneset = NULL, condition = "Adjust.P.Value", cutoff = 0.05, export = TRUE, path = NULL)|
| 3.5.2 Favorability Score                                           | `TimiFS`(resdata = NULL, condition = "Adjust.P.Value", cutoff = 0.05)|
| 3.5.3 Visualization: Bar plot of Favorability Score                | `TimiFSBar`(score = NULL, select = NULL)|

## Available Data
Here is a summary of available data in the package:
1. Cell type and marker annotation:

| Available Data                                             | Description                                                       |
| ---------------------------------------------------------- | ------------------------------------------------------------------| 
|    data(Immune_Marker_n1293)                            | #1293 immune cell markers                                         |
|    data(`CellType_Charoentong2017_Bindea2013_Xu2018_Immune`)| Immune cell types and markers from 3 publications             |
|    data(`CellType_Bindea2013_cancer`)                   | Markers of immune cell types and 1 cancer cell type               |
|    data(`CellType_Newman2015_LM22`)                     | Markers of immune cell types(LM22 signature) generated by CIBERSORT|
|    data(`CellType_Tirosh2016_melanoma_TME`)             | Metastatic melanoma tumor microenvironment cell types and markers |
|    data(`CellType_Zheng2021_Tcell`)                     | T cell subtypes and markers                                       |
|    data(`CellType_Zheng2021_Tumor`)                     | Tumor cell types and markers                                      |
|    data(`CellType_TNBC_aPDL1`)                     | Immune cell types and markers customized for TNBC (anti-PD1/PD-L1)                                      |

2. Example of Prognosis Module:

| Available Data                                             | Description                                                       |
| ---------------------------------------------------------- | ------------------------------------------------------------------| 
|   data(SKCM06info)                                     | Input: TCGA SKCM06(metastatic melanoma) clinical info                     |
|    data(SKCM06rna)                                      | Input: TCGA SKCM06(metastatic melanoma) transcriptome profile             |
| data(Bindea2013c_COX_MP_SKCM06)                      | Result: `TimiCOX` results with the `CellType_Bindea2013_cancer` annotation that reveals the association between each marker pair and favorable prognosis.   |
|  data(Bindea2013c_enrich)                             | Result: An example of `TimiEnrich` and `TimiPermFDR` results on cell interactions, analyzed following above          |

3. Example of Response Module:

| Available Data                                             | Description                                                       |
| ---------------------------------------------------------- | ------------------------------------------------------------------|
|   data(TNBCaPD1rna)                                   | Input: GSE194040 TNBC chemoimmunotherapy cohort microarray expression |
|   data(TNBCaPD1info)                                  | Input: GSE194040 TNBC chemoimmunotherapy cohort therapy response info |



## Example
Here is an example[01](example/example01_Bindea2013_Cancer.R) of how to use Prognosis Module. 
The useage of Response Module is similar to the Prognosis Module. Please refer to the example[07](example/example07_customized_TNBC_Response.R) for more details.
Other examples can be found in the [example](example/) folder. And the process to generate the cell markers can be found in the [inst/extdata](inst/extdata) folder.

Library the package
```R
library(TimiGP)
rm(list=ls())
```

### 1.   Preprocess of clinical information and transcriptome profile
`TimiCheckEvent` will filter out the low-quality clinical info(event and time-to-event).

`TimiPrePropress` will extract the transcription profile of selected markers and cohorts, then perform log transformation and gene-wise median normalization.
```R
rm(list=ls())
#1. Load SCKCM06 data ----
data("SKCM06info")
head(SKCM06info)
data("SKCM06rna")
#2. Load cell type and marker annotation ----
data("CellType_Bindea2013_cancer")
geneset <- CellType_Bindea2013_cancer
marker <- unique(geneset$Gene)
#3. Preprocess: TimiCheckEvent & TimiPrePropress ----
info <- TimiCheckEvent(SKCM06info)
rna <- TimiPrePropress(marker = marker,rna = SKCM06rna,cohort = rownames(info))
```
### 2 Gene Interaction
#### 2.1   Pairwise Comparison
In default, `TimiGenePair` will capture the logical relation of any two marker pairs, and generate a matrix of Marker Pair Score(MPS):

  - 1 or TRUE = the expression of gene A > that of gene B, 

  - 0 or FALSE = the expression of gene A < that of gene B.

Optional:  Capture the continuous relation of any two marker pairs by setting `cont = T`, and generate a matrix of Marker Pair Scores:

  -  the expression of gene A - that of gene B (See [example06_Bindea2013_Cancer_continuous_pair.R](example/example06_Bindea2013_Cancer_continuous_pair.R))
```R
#4. Generate marker pair score: TimiGenePair  ----
mps <- TimiGenePair(rna)
```
#### 2.2   Directed IMGP Selection
`TimiCOX` will perform univariate Cox regression that fits each marker pair as a variable. The result of Cox regression is returned as the first list. If a Pair A_B is associated with a poor prognosis(HR > 1), even if not significant, it will be changed to B_A and reverse its value in the matrix of the pair. The new matrix of Marker Pair Score(MPS) is returned as the second list.


This step takes about 5-10 min, which depends on the number of gene pairs.

```R
#5. Perform univariate Cox regression to find the association between marker pair and survival: TimiCOX ----
res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
mps <- res[[1]]
cox_res <- res[[2]]

# This step takes about 20-30 min, the result has been saved in data as examples
# Bindea2013c_COX_MP_SKCM06 <- cox_res
# save(Bindea2013c_COX_MP_SKCM06, file = "data/Bindea2013c_COX_MP_SKCM06.rda")
```
#### 2.3   Directed Gene Network
By setting "export = TRUE", `TimiGeneNetwork` generates three files that can be used to build a network in Cytoscape: 
 - network files: simple interaction file (gene_network.sif); 
 - node attributes (gene_node.txt); 
 - edge attributes (gene_edge.txt). 
The function also returns a list of the above files that can be modified in R.
```R
rm(list=ls())
# 6. Generate Directed Gene Network: TimiGeneNetwork  ----
data(Bindea2013c_COX_MP_SKCM06)
cox_res <- Bindea2013c_COX_MP_SKCM06
# You can use Cytoscape to visualize the network by setting "export = TRUE"
NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Bindea2013_Cancer",export =TRUE, path = "./")
```
### 3 Cell Interaction
`TimiCellPair` will generate marker pair annotations of cell interactions. For example, given any two different cell types, Cell A has markers a1 and a2 and Cell B has markers b1 and b2. Cell interaction A_B includes marker pairs: a1_b1, a1_b2, a2_b1, a2_b2 while cell interaction B_A includes marker pairs: b1_a1, b1_a2, b2_a1, b2_a2.
#### 3.1   Cell Interaction Annotation
```R
rm(list=ls())
# 7. Generate Cell Interaction Annotation: TimiCellPair ----
data(CellType_Bindea2013_cancer)
geneset <- CellType_Bindea2013_cancer
cell_pair <- TimiCellPair(geneset = geneset,core = 20)
```
#### 3.2   Prepare Enrichment Background
`TimiBG` will generate background gene pairs for enrichment analysis. The background should include both directions of each pair. For example, given a gene pair A_B, it will generate the background gene pairs A_B and B_A.
```R
# 8. Generate background: TimiBG ----
background <- TimiBG(marker.pair = row.names(cox_res))
```
#### 3.3   Enrichment Analysis
`TimiEnrich` will perform the analysis according to the query gene pairs. 
```R
# 9. Query: Select marker pairs A_B=1 significantly associated with good prognosis ----
data(Bindea2013c_COX_MP_SKCM06)
cox_res <- Bindea2013c_COX_MP_SKCM06
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]
# 10. Enrichment Analysis: TimiEnrich ----
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)
                  
# Optional: Permutation to generate FDR: TimiPermFFDR ----
res <- TimiPermFDR(resdata = res, geneset = geneset, gene = GP,
                   background = background, niter = 100, core = 20)
# This step takes some time. Please be patient

# This has been saved to data as an example
# Bindea2013c_enrich <- res
# save(Bindea2013c_enrich,file = "data/Bindea2013c_enrich.rda")
```
`TimiDotplot` can visualize the enrichment analysis results.
```R
# 11. Visualization: Dot plot of selected cell interaction: TimiDotplot-----
rm(list=ls())
data("Bindea2013c_enrich")
res <- Bindea2013c_enrich
p <- TimiDotplot(resdata = res,select = c(1:10))
p
```
![Fig3A](/assets/images/Fig3A.png)

`TimiCellChord` can visualize the cell interaction in Chord Diagram.
```R
# 12. Visualization: Chord Diagram of functional interaction: TimiCellChord----
rm(list=ls())
data("Bindea2013c_enrich")
res <- Bindea2013c_enrich
# Cell Chord Diagram
TimiCellChord(resdata = res,dataset = "Bindea2013_Cancer")
```
![Fig3B](/assets/images/Fig3B.png)

If you are interested in enriched marker pairs in specific cell interactions such as Cytotoxic cells_Cancer cells, please use `TimiGeneChord`.
```R
# Chord Diagram of marker pairs in select cell interaction
TimiGeneChord(resdata = res,select = 3)
```
![Fig3x](/assets/images/Fig3x.png)
#### 3.4  Cell Interaction Network
By setting "export = TRUE", `TimiCellNetwork` generates three files that can be used to build the network in Cytoscape: 
 - network files: simple interaction file (cell_network.sif); 
 - node attributes (cell_node.txt); 
 - edge attributes (cell_edge.txt). 
The function also returns a list of the above files that can be modified in R.

```R
# 13. Generate Directed Cell Network: TimiCellNetwork  ----
# You can use Cytoscape to visualize the network
rm(list=ls())
data("Bindea2013c_enrich")
res <- Bindea2013c_enrich
NET <- TimiCellNetwork(resdata = res,dataset = "Bindea2013_Cancer",export =TRUE, path = "./")
```
![Fig3C](/assets/images/Fig3C.png)
#### 3.5  Favorability Score
Based on the degree information of the cell interaction network, `TimiFS` calculates the favorability score of each cell type. In the network, 

 - the out-degree of cell A means how many times high A-over-other cell function is associated with favorable prognosis; 

 - the in-degree of cell A means how many times low A-over-other cell function is associated with favorable prognosis. 

The favorability score includes a favorable score and an unfavorable score. 
 - Favorable score: out-degree of the cell/sum of out-degree of all cell*100. 
 - Unfavorable score: in-degree of the cell/sum of in-degree of all cell*100.


```R
# 13. Calculate favorability score: TimiFS ----
# Visualization: TimiFSBar 
rm(list=ls())
data("Bindea2013c_enrich")
res <- Bindea2013c_enrich
# Calculate
score <- TimiFS(res)
head(score)
```
This result can be visualized by `TimiFSBar`.
```R
# Visualization
p <- TimiFSBar(score,select = c(1:5,(nrow(score)-2):nrow(score)))
p
```
![Fig3D](/assets/images/Fig3D.png)
