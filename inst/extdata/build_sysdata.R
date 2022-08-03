# sysdata are those markers
## setwd("./") # working directory is the TimiGP folder

rm(list = ls())

myinf1 <- "data/CellType_Charoentong2017_Bindea2013_Xu2018_Immune.rda"
myinf2 <- "data/CellType_Hughes2020_Tirosh2016_melanoma_TME.rda"
myinf3 <- "data/CellType_Zheng2021_Tcell.rda"
myinf4 <- "data/CellType_Bindea2013_cancer.rda"

load(myinf1)
load(myinf2)
load(myinf3)
load(myinf4)

Ann_Immune <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune
Ann_melanoma_TME <- CellType_Hughes2020_Tirosh2016_melanoma_TME
Ann_Tcell <-CellType_Zheng2021_Tcell
Ann_Bindea2013_Cancer <-  CellType_Bindea2013_cancer
# save(Ann_Immune, Ann_melanoma_TME,Ann_Tcell, Ann_Bindea2013_Cancer, file = "R/sysdata.rda")