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
| Step                                                               | Function       | Description                                                  |
| ------------------------------------------------------------------ | ---------- | ------------------------------------------------------------ |
| 1.1   Preprocess of Clinical Info                                  | TimiCheckEvent(SKCM06info) |  |
| 1.2   Preprocess of Transcriptome Profile                          | TimiPrePropress(marker = marker,rna = SKCM06rna,cohort = rownames(info))||
| 2.1   Pairwise Comparison                                          | TimiGenePair(rna)
| 2.2   Directed IRGP Selection                                      | TimiCOX(mps = mps,info = info,p.adj = "BH")
| 2.3   Directed Gene Network                                        | TimiGeneNetwork(resdata = cox_res,dataset = "Galon2013_Cancer",path = "./")
| 3.1   Cell Pair Annotation                                         | TimiCellPair(geneset = geneset,core = 20)
| 3.2   Prepare Enrichment Background                                | TimiBG(marker.pair = row.names(cox_res))
| 3.3.1 Enrichment Analysis                                          | TimiEnrich(gene = GP, background = background, geneset = cell_pair, p.adj = "BH",core=20)
| 3.3.2 Visualization: Dotplot of Enrichment                         | TimiDotplot(resdata = res,select = c(1:10))
| 3.4.1 Visualization: Chord Diagram of Cell Interaction             | TimiCellChord(resdata = res,dataset = "Galon2013_Cancer")
| 3.4.2 Visualization: Chord Diagram of Selected Gene Interaction    | TimiGeneChord(resdata = res,select = 1)
| 3.5   Cell Interaction Network                                     | TimiCellNetwork(resdata = res,dataset = "Galon2013_Cancer",path = "./")
| 3.6   Favorability Score                                           | TimiFS(res)
| 3.7   Visualization: Bar plot of Favorability Score                | TimiFSBar(score)


