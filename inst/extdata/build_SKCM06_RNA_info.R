# The Cancer Genome Atlas (TCGA) RNA-seq dataset for skin cutaneous melanoma (SKCM) 
# was downloaded from Firehose (https://gdac.broadinstitute.org/). 
# This dataset consisted of RSEM normalized gene expression data for 20,501 genes from 473 skin cutaneous melanoma samples. 
# Of the SKCM samples, 103 are from primary tumors and 368 are from metastatic tumor samples. 
# We used the metastatic tumor samples --- SKCM06

## SKCM06 clinical info#####
myinf1 <- "~/TCGA_SKCM/SKCM_RNAseqv2_ALL_Symbol.rda"
myinf2 <- "~/TCGA_SKCM/SKCM_Clincial_info.txt"

myoutf1 <- "./data/SKCM06rna.rda"
myoutf2 <- "data/SKCM06info.rda"
## load SKCM06 RNA
load(myinf1)
rna <- mydata
xx <- as.numeric(substr(colnames(rna), 14, 15))
se <- which(xx==6)		## metastatic
length(se) # 368 samples
rna <- rna[,se]
dim(rna) # 20501   368
colnames(rna) <- substr(colnames(rna), 1, 12)

SKCM06rna <- rna
# save(SKCM06rna,file = myoutf1) 

## load clinical info
info <- read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
head(info)
colnames(info)

se <- c("vital_status", "gender", "stage_event.pathologic_stage", "days_to_death", "days_to_last_followup",  "age_at_initial_pathologic_diagnosis", "race_list.race", "tumor_tissue_site")
info <- info[,se]

xx <- info[, "vital_status"]
event <- ifelse(xx=="dead", 1, 0)
event[is.na(xx)] <- NA
xx <- info[, "days_to_death"] 
time <- ifelse(!is.na(xx), xx, info[, "days_to_last_followup"])
time <- as.numeric(time)
SKCM06info <- data.frame( event, time)
rownames(SKCM06info ) <- rownames(info)

comSam <- intersect(row.names(SKCM06info), colnames(rna))


SKCM06info <- SKCM06info[comSam,]
# save(SKCM06info ,file = myoutf2)
