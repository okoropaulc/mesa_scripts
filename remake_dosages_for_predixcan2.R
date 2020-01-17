#check the en ALL model summaries, and see why using it to impute gives so many NA's
#ryan has the ALL model summaries without header, so read in summaries with header and copy the col
library(data.table)
# en_cau <- fread(file="Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T, nrows=3)
# en_all <- fread(file="Z:/no_header_MESA.ALL.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=F)
# 
# header <- colnames(en_cau)
# 
# names(en_all) <- header
# 
# rf_all <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
# 
# en_all <- subset(en_all, cv_R2_avg > 0.01)
# rf_all <- subset(rf_all, CV_R2 > 0.01)
# 

#check the initial dosage
# cau_dos <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/chr22txt.gz", header=T, nrows=5)
# cau_dos[1:5,1:10]
# mycau_dos <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr22_all_sample.txt", header=T)
# mycau_dos[1:5,1:10]
# impcau <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr22.txt", header=T)
# impcau[1:5,1:10]
# 
# imp37 <- as.data.frame(impcau$id)
# colnames(imp37) <- "b37"
# imp37$b37 <- as.character(imp37$b37)
# 

#pass arguments to the script
args <- commandArgs(trailingOnly = T)
chrom <- as.character(args[1])

"%&%" <- function(a,b) paste(a,b, sep = "")

library(data.table)
print(chrom)

#check the weight file, and use it to change to dosages to have same rsid
allwei <- fread(file="/home/pokoro/MESA.ALL.WG.PC3.PF10.unpruned.rsid.hg19.weights.txt", header=F)
#allwei[1:5,1:6]

allwei37 <- as.data.frame(allwei$V3)
colnames(allwei37) <- "b37"
allwei37$b37 <- as.character(allwei37$b37)

library(dplyr)

#b37join <- inner_join(imp37, allwei37, by = c("b37" = "b37")) #7147 b37 overlaps. therefore ryan cpos liftover is wrong

# check rsid overlap
#allwei$V3 <- as.character(allwei$V3)
#rsid_join <- inner_join(allwei, b37join, by = c("V3" = "b37"))

library(tidyverse)
#rsid_join <- drop_na(rsid_join)
#rsid_join <- unique(rsid_join) #7147 rsids overlap

#therefore use rsid in the dosage (that is, remake the dosage), and use the ALL_unfiltered_rsid_hg19 build to predict expression

cau_dos <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan/chr" %&% chrom %&% "_cau.txt.gz", header=F)
#cau_dos[1:5,1:10]

#cau_dos1 <- cau_dos[1:5000,1:10]
cau_dos$V1 <- as.character(cau_dos$V1)
cau_dos$V2 <- as.character(cau_dos$V2)
cau_dos$V3 <- as.character(cau_dos$V3)
cau_dos$V4 <- as.character(cau_dos$V4)
cau_dos$V5 <- as.character(cau_dos$V5)

for (j in 1:nrow(cau_dos)){
  cau_dos$V2[j] <- as.character(cau_dos$V1[j] %&% "_" %&% cau_dos$V3[j] %&% "_" %&% cau_dos$V4[j] %&% "_" %&% cau_dos$V5[j] %&% "_b37")
}

rsid2 <- inner_join(allwei, cau_dos, by = c("V3"="V2"))
rsid2 <- rsid2[,c(2,7:length(rsid2))] #take out rsid and the other dosage columns ie chr pos ref alt maf id1....idn
rsid2 <- rsid2[,c(2,1,3:length(rsid2))]
rsid2 <- unique(rsid2) #remove duplicates
rsid2[1:5,1:10]

#write out without column names. chromosome rsid position allele1 allele2 MAF id1 ..... idn. predixcan does not need colnames
fwrite(rsid2, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan/correct_dosages/chr" %&% chrom %&% "_cau.txt", quote=F, row.names=F, sep ="\t", col.names=F)
