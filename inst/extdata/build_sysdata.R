# sysdata are those markers
## setwd("./") # working directory is the TimiGP folder

rm(list = ls())
data("CellType_Charoentong2017_Galon2013_TIP_Immune")
Ann_Immune <- CellType_Charoentong2017_Galon2013_TIP_Immune
data("CellType_Alex2020_Levi2019_TME")
Ann_TME <- CellType_Alex2020_Levi2019_TME
data("CellType_Zhang2021S2_Tcell")
Ann_Tcell <-CellType_Zhang2021S2_Tcell
data("CellType_Galon2013_cancer")
Ann_Galon2013_Cancer <-  CellType_Galon2013_cancer
# save(Ann_Immune, Ann_TME,Ann_Tcell, Ann_Galon2013_Cancer, file = "R/sysdata.rda")