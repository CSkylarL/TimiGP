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


