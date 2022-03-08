#' Cell Type in melanoma microenvironment generated from Alex2020 and Levi2019 with cancer cells
#'
#' A dataset containing the cell type and markers annotated by Alex2020 and Levi2019
#' 
#' @docType data
#'
#' @usage data(CellType_Alex2020_Levi2019_TME)
#' 
#' @keywords datasets
#' 
#' @format A data frame with 663 rows and 3 variables:
#' \describe{
#'   \item{CellType}{Immune cell type, skin cell type and cancer cell type}
#'   \item{Gene}{Marker gene of the cell type}
#'   \item{Dataset}{The source of the annotation}
#' }
#' 
#' @references Hughes, T. K., Wadsworth II, M. H., Gierahn, T. M., Do, T., Weiss, D., Andrade, P. R., ... & Shalek, A. K. (2020). Second-strand synthesis-based massively parallel scRNA-seq reveals cellular states and molecular features of human inflammatory skin pathologies. Immunity, 53(4), 878-894.
#' @references Tirosh, I., Izar, B., Prakadan, S. M., Wadsworth, M. H., Treacy, D., Trombetta, J. J., ... & Garraway, L. A. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science, 352(6282), 189-196.
#'
#' @source Modification can be found in inst/extdata/build_CellType_Alex2020_Levi2019_TME.R
"CellType_Alex2020_Levi2019_TME"







#' Cell Type generated from Galon2013 with cancer cells
#'
#' A dataset containing the cell type and markers annotated by Galon2013
#' 
#' @docType data
#'
#' @usage data(CellType_Galon2013_cancer)
#' 
#' @keywords datasets
#' 
#' @format A data frame with 579 rows and 3 variables:
#' \describe{
#'   \item{CellType}{Immune cell type and cancer cell type}
#'   \item{Gene}{Marker gene of the cell type}
#'   \item{Dataset}{The source of the annotation}
#' }
#' 
#' @references Bindea, G., Mlecnik, B., Tosolini, M., Kirilovsky, A., Waldner, M., 
#' Obenauf, A. C., ... & Galon, J. (2013). Spatiotemporal dynamics of intratumoral 
#' immune cells reveal the immune landscape in human cancer. Immunity, 39(4), 782-795.
#'
#' @source Modification can be found in inst/extdata/build_CellType_Galon2013_cancer.R
"CellType_Galon2013_cancer"








#' Immune Cell Type generated from Charoentong2017, Galon2013, TIP
#'
#' A dataset containing the cell type and markers annotated by Charoentong2017, Galon2013, TIP
#' 
#' @docType data
#'
#' @usage data(CellType_Charoentong2017_Galon2013_TIP_Immune)
#' 
#' @keywords datasets
#' 
#' @format A data frame with 1409 rows and 3 variables:
#' \describe{
#'   \item{CellType}{Immune cell type}
#'   \item{Gene}{Marker gene of the cell type}
#'   \item{Dataset}{The source of the annotation}
#' }
#' 
#' @references Charoentong, P., Finotello, F., Angelova, M., Mayer, C., Efremova, M., Rieder, D., ... & Trajanoski, Z. (2017). Pan-cancer immunogenomic analyses reveal genotype-immunophenotype relationships and predictors of response to checkpoint blockade. Cell reports, 18(1), 248-262.
#' @references Bindea, G., Mlecnik, B., Tosolini, M., Kirilovsky, A., Waldner, M., Obenauf, A. C., ... & Galon, J. (2013). Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity, 39(4), 782-795.
#' @references Xu, L., Deng, C., Pang, B., Zhang, X., Liu, W., Liao, G., ... & Li, X. (2018). TIP: a web server for resolving tumor immunophenotype profiling. Cancer research, 78(23), 6575-6580.
#'
#' @source Modification can be found in inst/extdata/build_CellType_Charoentong2017_Galon2013_TIP_Immune.R
"CellType_Charoentong2017_Galon2013_TIP_Immune"









#' T Cell Subtypes generated from Zhang2021S2
#'
#' A dataset containing the cell type and markers annotated by Zhang2021S2
#' 
#' @docType data
#'
#' @usage data(CellType_Zhang2021S2_Tcell)
#' 
#' @keywords datasets
#' 
#' @format A data frame with 926 rows and 3 variables:
#' \describe{
#'   \item{CellType}{T cell subtype}
#'   \item{Gene}{Marker gene of the cell type}
#'   \item{Dataset}{The source of the annotation}
#' }
#' 
#' @references Zheng, L., Qin, S., Si, W., Wang, A., Xing, B., Gao, R., ... & Zhang, Z. (2021). Pan-cancer single-cell landscape of tumor-infiltrating T cells. Science, 374(6574), abe6474. 
#'
#' @source Modification can be found in inst/extdata/build_CellType_Zhang2021S2_Tcell.R
"CellType_Zhang2021S2_Tcell"






#' Immune markers generated from Charoentong2017, Galon2013, TIP and classic markers
#'
#' A dataset containing the immune markers annotated by Charoentong2017, Galon2013, TIP 
#' and additional classic immune inhibitors, stimulators and checkpoints.
#' 
#' @docType data
#'
#' @usage data(Immune_Marker_n1326)
#' 
#' @keywords datasets
#' 
#' @format A vector of 1326 immune markers
#' 
#' @references Charoentong, P., Finotello, F., Angelova, M., Mayer, C., Efremova, M., Rieder, D., ... & Trajanoski, Z. (2017). Pan-cancer immunogenomic analyses reveal genotype-immunophenotype relationships and predictors of response to checkpoint blockade. Cell reports, 18(1), 248-262.
#' @references Bindea, G., Mlecnik, B., Tosolini, M., Kirilovsky, A., Waldner, M., Obenauf, A. C., ... & Galon, J. (2013). Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity, 39(4), 782-795.
#' @references Xu, L., Deng, C., Pang, B., Zhang, X., Liu, W., Liao, G., ... & Li, X. (2018). TIP: a web server for resolving tumor immunophenotype profiling. Cancer research, 78(23), 6575-6580.
#'
#' @source Modification can be found in inst/extdata/build_Immune_Marker_n1326.R
"Immune_Marker_n1326"








#' TCGA metastatic melanoma(SKCM06) RNA expression
#'
#' A dataset containing 20,501 genes from 368  metastatic melanoma samples.
#' 
#' @docType data
#'
#' @usage data(SKCM06rna)
#' 
#' @keywords datasets
#' 
#' @format A data frame with 20501 rows and 368 variables:
#' \describe{
#'   \item{row}{Gene}
#'   \item{column}{Sample(patient) ID}
#' }
#' 
#' 
#' @references \url{"https://gdac.broadinstitute.org/"}
#' 
#' @source Modification can be found in inst/extdata/build_SKCM06_RNA_info.R
"SKCM06rna"













#' TCGA metastatic melanoma(SKCM06) clinical infomation
#'
#' A dataset containing 20,501 genes from 368  metastatic melanoma samples.
#' 
#' @docType data
#'
#' @usage data(SKCM06info)
#' 
#' @keywords datasets
#' 
#' @format A data frame with 368 rows and 2 variables:
#' \describe{
#'   \item{row}{Patient ID}
#'   \item{column}{Clinical information, including event(death) and time(day) to event}
#' }
#' 
#' 
#' @references \url{"https://gdac.broadinstitute.org/"}
#' 
#' @source Modification can be found in inst/extdata/build_SKCM06_RNA_info.R
"SKCM06info"




#' COX regression Results from function TimiCOX with cell type meaker annotated by Galon2013_cancer
#'
#' An intermediate result generated from function TimiCOX
#' that reveals the association between each marker pairs and favorable prognosis.
#' 
#' @docType data
#'
#' @usage data(Galon2013c_COX_MP_SKCM06)
#' 
#' @keywords intermediate result
#' 
#' @format A data frame with 127765 rows and 3 variables:
#' \describe{
#'   \item{Row name}{Marker pair}
#'   \item{HR}{Hazard.Ratio}
#'   \item{PV}{P-Value}
#'   \item{QV}{Adjust P-value}
#' }
#' 
#' 
#' @source intermediate result generated from function TimiCOX
"Galon2013c_COX_MP_SKCM06"









#' COX regression Results from function TimiCOX with cell type meaker annotated by Charoentong2017_Galon2013_TIP_Immune
#'
#' An intermediate result generated from function TimiCOX
#' that reveals the association between each marker pairs and favorable prognosis.
#' 
#' @docType data
#'
#' @usage data(Immune3_COX_MP_SKCM06)
#' 
#' @keywords intermediate result
#' 
#' @format A data frame with 739936 rows and 3 variables:
#' \describe{
#'   \item{Row name}{Marker pair}
#'   \item{HR}{Hazard.Ratio}
#'   \item{PV}{P-Value}
#'   \item{QV}{Adjust P-value}
#' }
#' 
#' 
#' @source intermediate result generated from function TimiCOX
"Immune3_COX_MP_SKCM06"
















#' COX regression Results from function TimiCOX with cell type meaker annotated by Alex2020_Levi2019_TME
#'
#' An intermediate result generated from function TimiCOX
#' that reveals the association between each marker pairs and favorable prognosis.
#' 
#' @docType data
#'
#' @usage data(TME_COX_MP_SKCM06)
#' 
#' @keywords intermediate result
#' 
#' @format A data frame with 153735 rows and 3 variables:
#' \describe{
#'   \item{Row name}{Marker pair}
#'   \item{HR}{Hazard.Ratio}
#'   \item{PV}{P-Value}
#'   \item{QV}{Adjust P-value}
#' }
#' 
#' 
#' @source intermediate result generated from function TimiCOX
"TME_COX_MP_SKCM06"















#' COX regression Results from function TimiCOX with cell type meaker annotated by Zhang2021S2_Tcell
#'
#' An intermediate result generated from function TimiCOX
#' that reveals the association between each marker pairs and favorable prognosis.
#' 
#' @docType data
#'
#' @usage data(Tcell_COX_MP_SKCM06)
#' 
#' @keywords intermediate result
#' 
#' @format A data frame with 153735 rows and 3 variables:
#' \describe{
#'   \item{Row name}{Marker pair}
#'   \item{HR}{Hazard.Ratio}
#'   \item{PV}{P-Value}
#'   \item{QV}{Adjust P-value}
#' }
#' 
#' 
#' @source intermediate result generated from function TimiCOX
"Tcell_COX_MP_SKCM06"


































#' Cell Interaction Enrichment from function TimiEnrich 
#' with cell type meaker annotated by Galon2013_cancer
#'
#' An example result of cell interaction enrichment generated from function TimiEnrich
#' 
#' @docType data
#'
#' @usage data(Galon2013c_enrich)
#' 
#' @keywords example result
#' 
#' 
#' @format a data frame of enrichment result:
#' \describe{
#'   \item{Index}{A numeric vector used to specify which cell pair would be used in visualization}
#'   \item{Rank}{A rank of cell interaction according to Adjust. P-value}
#'   \item{Cell.Pair}{Cell Interaction, e.g. Cell A_Cell B}
#'   \item{Favorable.Cell.Type}{In the cell pair A_B, the first cell A}
#'   \item{Unfavorable.Cell.Type}{In the cell pair A_B, the second cell B}
#'   \item{No.Total.Marker.Pair}{The total number of marker pairs in the cell pair annotation}
#'   \item{No.Shared.Marker.Pair}{The number of marker pairs shared by query pairs and annotation}
#'   \item{Enrichment.Ratio}{Enrichment Ratio}
#'   \item{P.Value}{P. value}
#'   \item{Adjust.P.Value}{Adjust.P.Value}
#'   \item{Shared.Marker.Pair}{The marker pairs shared by query pairs and annotation}
#'   \item{Total.Marker.Pair}{The marker pairs in the annotation}
#' }
#' 
#' 
#' @source example result generated from function TimiEnrich(example/example01_Galon2013_Cancer.R)
"Galon2013c_enrich"
