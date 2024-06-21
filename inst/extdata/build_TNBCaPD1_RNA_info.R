# Data source ##################################################################
# The microarray dataset GSE194040 for Triple-negative breast cancer (TNBC)
# was downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194040&format=file&file=GSE194040%5FISPY2ResID%5FAgilentGeneExp%5F990%5FFrshFrzn%5FmeanCol%5FgeneLevel%5Fn988%2Etxt%2Egz). 
# The file was decompressed and saved to `myinf1`

# The corresponding clinical info
# was downloaded from GEO (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194040/soft/GSE194040_family.soft.gz)
# The file was decompressed and saved to 
# "./Breast_Cancer/I-SPY_GSE194040/raw/GSE194040_family.soft"

# The  clinical info was preprocessed with below shell script
# And saved to `myinf2`
# sh Extract_clinical_info_from_familysort.sh -->
#!/bin/sh
# cat *family.soft |grep ! |more
# cat *family.soft |grep SAMPLE|cut -d "=" -f 2|sed s#"^ "##g > sample_tmp.txt
# cat *family.soft |grep "patient id:"|cut -d ":" -f 2|sed s#"^ "##g > title_tmp.txt
# cat *family.soft |grep Sample_source_name_ch1|cut -d "=" -f 2|sed s#"^ "##g > tissue_tmp.txt
# cat *family.soft |grep "her2:"|cut -d ":" -f 2|sed s#"^ "##g > her2_tmp.txt
# cat *family.soft |grep "hr:"|cut -d ":" -f 2|sed s#"^ "##g > hr_tmp.txt
# cat *family.soft |grep "pcr:"|cut -d ":" -f 2|sed s#"^ "##g > pcr_tmp.txt
# cat *family.soft |grep "mp:"|cut -d ":" -f 2|sed s#"^ "##g > mp_tmp.txt
# cat *family.soft |grep "arm:"|cut -d ":" -f 2|sed s#"^ "##g > arm_tmp.txt
# ls *tmp.txt|cut -d "_" -f 1| sed ':a;N;$!ba;s#\n#\t#g' > Clinical_info.txt
# paste -d "\t" *tmp.txt >> Clinical_info.txt
# mv Clinical_info.txt ../
# rm *tmp.txt

# clear workspace and load libraries 
rm(list = ls())
library(dplyr)
library(stringr)
library(data.table)
library(tibble)
# Filter data ################################################################
myinf1 <- "./Breast_Cancer/I-SPY_GSE194040/raw/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt"
myinf2 <- "./Breast_Cancer/I-SPY_GSE194040/Clinical_info.txt"

myoutf1 <- "./data/TNBCaPD1rna.rda"
myoutf2 <- "data/TNBCaPD1info.rda"

# rna-----------------------------------------------------------------------
rna <- read.table(myinf1,sep = "\t",header = T,row.names = 1)

# info----------------------------------------------------------------------
info <-  read.table(myinf2,sep = "\t",header = T)

table(info$pcr, exclude = F) # NA 1 drop

# Preprocess the clinical information 
info <- info %>% 
  mutate(rownames = paste0("X",title)) %>%
  filter(rownames %in% colnames(rna)) %>%
  filter(!is.na(pcr)) %>% # remove unknown treatment response
  filter(arm == "Paclitaxel + Pembrolizumab") %>% # select the treatment
  filter(her2 == 0 & hr == 0) %>%  # filter TNBC
  mutate(Response =  pcr) %>%
  column_to_rownames(var = "rownames") %>%
  select( "Response" )

# filter the rna samples 
rna <- rna[ rownames(info)]

all(rownames(info) == colnames(rna)) # T

table(info$Response,  exclude = F)
# 0  1 
# 13 13 

# Save data ################################################################
TNBCaPD1rna <- rna
TNBCaPD1info <- info
# save(TNBCaPD1rna ,file = myoutf1)
# save(TNBCaPD1info ,file = myoutf2)
