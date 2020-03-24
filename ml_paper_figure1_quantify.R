hdl <- fread(file="Z:/data/twas_mesa/hdl_notrain.txt", header = T)

library(dplyr)

afa <- subset(afa_sam, V2 %in% hdl$sidno)
his <- subset(his_sam, V2 %in% hdl$sidno)
cau <- subset(cau_sam, V2 %in% hdl$sidno)

###########               ALL

en_all <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", header=T)#it has no header
en_all <- en_all[,c(1,2,10)]
library(tidyverse)
en_all <- drop_na(en_all) #remove NA
en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)
#en_all <- subset(en_all, cv_R2_avg > -0.5)
en_all$gene_id <- as.character(en_all$gene_id)


rf_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
rf_all$CV_R2 <- as.numeric(rf_all$CV_R2)
#rf_all <- subset(rf_all, CV_R2 > -0.5)
rf_all$Gene_ID <- as.character(rf_all$Gene_ID)


svr_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
svr_all$CV_R2 <- as.numeric(svr_all$CV_R2)
#svr_all <- subset(svr_all, CV_R2 > -0.5)
svr_all$Gene_ID <- as.character(svr_all$Gene_ID)


knn_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)
knn_all$CV_R2 <- as.numeric(knn_all$CV_R2)
#knn_all <- subset(knn_all, CV_R2 > -0.5)
knn_all$Gene_ID <- as.character(knn_all$Gene_ID)


library(dplyr)
en_rf <- inner_join(en_all, rf_all, by = c("gene_id"="Gene_ID"))
en_rf <- en_rf[,c(3,5)]
names(en_rf) <- c("en", "rf")


en_svr <- inner_join(en_all, svr_all, by = c("gene_id"="Gene_ID"))
en_svr <- en_svr[,c(3,5)]
names(en_svr) <- c("en", "svr")


en_knn <- inner_join(en_all, knn_all, by = c("gene_id"="Gene_ID"))
en_knn <- en_knn[,c(3,5)]
names(en_knn) <- c("en","knn")


#check number of genes EN performed better compared to the ML

knnvsen <- subset(en_knn, en_knn$knn > en_knn$en)
envsknn <- subset(en_knn, en_knn$en > en_knn$knn)
mean(envsknn$en) - mean(envsknn$knn)

rfvsen <- subset(en_rf, en_rf$rf > en_rf$en)
envsrf <- subset(en_rf, en_rf$en > en_rf$rf)
mean(envsrf$en) - mean(envsrf$rf)

svrvsen <- subset(en_svr, en_svr$svr > en_svr$en)
envssvr <- subset(en_svr, en_svr$en > en_svr$svr)
mean(envssvr$en) - mean(envssvr$svr)

#filter out R2 < 0.01
en_rf_0.1 <- subset(en_rf, en_rf$en > 0.01 & en_rf$rf > 0.01)
en_svr_0.1 <- subset(en_svr, en_svr$en > 0.01 & en_svr$svr > 0.01)
en_knn_0.1 <- subset(en_knn, en_knn$en > 0.01 & en_knn$knn > 0.01)


knnvsen_0.01 <- subset(en_knn_0.1, en_knn_0.1$knn > en_knn_0.1$en)
envsknn_0.01 <- subset(en_knn_0.1, en_knn_0.1$en > en_knn_0.1$knn)

rfvsen_0.01 <- subset(en_rf_0.1, en_rf_0.1$rf > en_rf_0.1$en)
envsrf_0.01 <- subset(en_rf_0.1, en_rf_0.1$en > en_rf_0.1$rf)

svrvsen_0.01 <- subset(en_svr_0.1, en_svr_0.1$svr > en_svr_0.1$en)
envssvr_0.01 <- subset(en_svr_0.1, en_svr_0.1$en > en_svr_0.1$svr)
