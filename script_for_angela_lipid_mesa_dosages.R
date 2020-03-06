#script to check the dosage and sample of angela's lipid paper

# pheno <- read.csv(file="Z:/data/mesa_models/mesa_pheno/Exam1Main.csv", header=T)
# chr1 <- read.table(file="Z:/data/mesa_models/mesa_pheno/angela_dosages/dosages1.txt.gz", header=F, nrows=3)
# sidno <- read.table(file="Z:/data/mesa_models/mesa_pheno/angela_dosages/samples.txt", header=F)
# 
# #keep only the sidno and hdl
# pheno <- pheno[,c(2,6)]
# 
# #check and remove any NA
# library(tidyverse)
# 
# #keep only pheno id that is in the sample id
# sidno$V2 <- as.character(sidno$V2)
# pheno$sidno <- as.character(pheno$sidno)
# 
# library(dplyr)
# 
# sample_pheno <- inner_join(sidno, pheno, by = c("V2"="sidno"))
# 
# sample_pheno <- sample_pheno[,c(2,8)]
# sample_pheno <- drop_na(sample_pheno)

#Use the mesa dosages from Lauren
library(data.table)
library(tidyverse)
library(dplyr)
pheno <- read.csv(file="Z:/data/mesa_models/mesa_pheno/Exam1Main.csv", header=T)
pheno$sidno<- as.character(pheno$sidno)

#Read in AFA
#afa_dos <- read.table(file="Z:/data/lauren_mesa/afa_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=3)
afa_sam <- fread(file="Z:/data/lauren_mesa/afa_dosages/samples.txt", header=F)

#Read in CAU
#cau_dos <- read.table(file="Z:/data/lauren_mesa/cau_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=3)
cau_sam <- fread(file="Z:/data/lauren_mesa/cau_dosages/samples.txt", header=F)

#Read in CAU
#his_dos <- read.table(file="Z:/data/lauren_mesa/his_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=3)
his_sam <- fread(file="Z:/data/lauren_mesa/his_dosages/samples.txt", header=F)

#phenotypes
trig_pheno <- pheno[,c(2,4)]
trig_pheno <- drop_na(trig_pheno)

ldl_pheno <- pheno[,c(2,5)]
ldl_pheno <- drop_na(ldl_pheno)

hdl_pheno <- pheno[,c(2,6)]
hdl_pheno <- drop_na(hdl_pheno)

chol_pheno <- pheno[,c(2,7)]
chol_pheno <- drop_na(chol_pheno)

#combine afa and cau samples
mesa_sam <- rbind(afa_sam, cau_sam, his_sam)
names(mesa_sam) <- c("FID","IID")
mesa_sam$FID <- NULL
mesa_sam$IID <- as.character(mesa_sam$IID)

#TWAS sample size per phenotype. That is, samples in the pheno and in the combined MESA
trig_mesa <- inner_join(trig_pheno, mesa_sam, by = c("sidno"="IID"))
ldl_mesa <- inner_join(ldl_pheno, mesa_sam, by = c("sidno"="IID"))
hdl_mesa <- inner_join(hdl_pheno, mesa_sam, by = c("sidno"="IID"))
chol_mesa <- inner_join(chol_pheno, mesa_sam, by = c("sidno"="IID"))
write.table(trig_mesa, file="Z:/data/twas_mesa/trig.txt",quote=F, row.names=F, sep="\t")
write.table(ldl_mesa, file="Z:/data/twas_mesa/ldl.txt",quote=F, row.names=F, sep="\t")
write.table(hdl_mesa, file="Z:/data/twas_mesa/hdl.txt",quote=F, row.names=F, sep="\t")
write.table(chol_mesa, file="Z:/data/twas_mesa/chol.txt",quote=F, row.names=F, sep="\t")


#remove samples used in the training from the phenos
# train_all_dos <- read.table(file="Z:/data/mesa_models/all/whole_genotypes/ALL.chr22.genotype.txt.gz",header=T,nrows=3,stringsAsFactors=F)
# train_all_sam <- colnames(train_all_dos)[2:length(train_all_dos)]
# for (i in 1:length(train_all_sam)){
#   train_all_sam[i] <- strsplit(train_all_sam[i], "X")[[1]][2]
# } #just to remove the X in from of the sample ID
# train_all_sam <- data.frame(FID=rep(0, length(train_all_sam)), IID=train_all_sam)
# #save the sample to file
# write.table(train_all_sam, file="Z:/data/mesa_models/all/whole_genotypes/samples.txt", col.names=F, row.names=F, sep="\t", quote=F)

train_all_sam <- fread(file="Z:/data/mesa_models/all/whole_genotypes/samples.txt", header=F)
names(train_all_sam) <- c("FID","IID")
train_all_sam$FID <- NULL
train_all_sam$IID <- as.character(train_all_sam$IID)

#Now remove the samples used for training from the phenos
trig_mesa_notrain <- anti_join(trig_mesa, train_all_sam, by = c("sidno"="IID"))
ldl_mesa_notrain <- anti_join(ldl_mesa, train_all_sam, by = c("sidno"="IID"))
hdl_mesa_notrain <- anti_join(hdl_mesa, train_all_sam, by = c("sidno"="IID"))
chol_mesa_notrain <- anti_join(chol_mesa, train_all_sam, by = c("sidno"="IID"))
write.table(trig_mesa_notrain, file="Z:/data/twas_mesa/trig_notrain.txt",quote=F, row.names=F, sep="\t")
write.table(ldl_mesa_notrain, file="Z:/data/twas_mesa/ldl_notrain.txt",quote=F, row.names=F, sep="\t")
write.table(hdl_mesa_notrain, file="Z:/data/twas_mesa/hdl_notrain.txt",quote=F, row.names=F, sep="\t")
write.table(chol_mesa_notrain, file="Z:/data/twas_mesa/chol_notrain.txt",quote=F, row.names=F, sep="\t")


#mesa samples for platelet count
thrombomodulin <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)
thrombomodulin$sidno <- as.character(thrombomodulin$sidno)
thrombomodulin_rankstm1 <- thrombomodulin[,c(1,4)]
thrombomodulin_rankstm1 <- drop_na(thrombomodulin_rankstm1)
thrombomodulin_rankplt5 <- thrombomodulin[,c(1,5)]
thrombomodulin_rankplt5 <- drop_na(thrombomodulin_rankplt5)

#TWAS sample size
rankstm1 <- inner_join(thrombomodulin_rankstm1, mesa_sam, by = c("sidno"="IID"))
rankplt5 <- inner_join(thrombomodulin_rankplt5, mesa_sam, by = c("sidno"="IID"))
write.table(rankplt5, file="Z:/data/twas_mesa/rankplt5.txt",quote=F, row.names=F, sep="\t")
write.table(rankstm1, file="Z:/data/twas_mesa/rankstm1.txt",quote=F, row.names=F, sep="\t")

#Remove train sample from the phenos
rankstm1_notrain <- anti_join(rankstm1, train_all_sam, by = c("sidno"="IID"))
rankplt5_notrain <- anti_join(rankplt5, train_all_sam, by = c("sidno"="IID"))
write.table(rankplt5_notrain, file="Z:/data/twas_mesa/rankplt5_notrain.txt",quote=F, row.names=F, sep="\t")
write.table(rankstm1_notrain, file="Z:/data/twas_mesa/rankstm1_notrain.txt",quote=F, row.names=F, sep="\t")

