#Make Figure 2 and 3
#Figure 2

#filter the mesa models with R2 and filter the spearman in MESA 2 METS
library(dplyr)
#####AFHI
en_afhi <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFHI_model_summaries.txt", header=T)#it has no header
en_afhi$gene_id <- as.character(en_afhi$gene_id)
en_afhi <- subset(en_afhi, en_afhi$cv_R2_avg > 0.01)
en_afhi <- en_afhi[,c(1,10)]
for (i in 1:length(en_afhi$gene_id)){
  en_afhi$gene_id[i] <- gsub('\\.[0-9]+','',en_afhi$gene_id[i])
} #just to remove the decimal places in the gene_id

en_afhi_2_mets <- read.table(file = "Z:/data/mesa_models/mets_dosages/new_predixcan_AFHI2METS_cor", header = T)
en_afhi_2_mets$gene <- as.character(en_afhi_2_mets$gene)
filt_en_afhi_2_mets <- inner_join(en_afhi, en_afhi_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01


#####ALL
en_all <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", header=T)#it has no header
en_all$gene_id <- as.character(en_all$gene_id)
en_all <- subset(en_all, en_all$cv_R2_avg > 0.01)
en_all <- en_all[,c(1,10)]
for (i in 1:length(en_all$gene_id)){
  en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
} #just to remove the decimal places in the gene_id

en_all_2_mets <- read.table(file = "Z:/data/mesa_models/mets_dosages/new_predixcan_ALL2METS_cor", header = T)
en_all_2_mets$gene <- as.character(en_all_2_mets$gene)
filt_en_all_2_mets <- inner_join(en_all, en_all_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01

#CAU
en_cau <- read.table(file = "Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
en_cau$gene_id <- as.character(en_cau$gene_id)
en_cau <- subset(en_cau, en_cau$cv_R2_avg > 0.01)
en_cau <- en_cau[,c(1,10)]
for (i in 1:length(en_cau$gene_id)){
  en_cau$gene_id[i] <- gsub('\\.[0-9]+','',en_cau$gene_id[i])
} #just to remove the decimal places in the gene_id

en_cau_2_mets <- read.table(file = "Z:/data/mesa_models/mets_dosages/new_predixcan_CAU2METS_cor", header = T)
en_cau_2_mets$gene <- as.character(en_cau_2_mets$gene)
filt_en_cau_2_mets <- inner_join(en_cau, en_cau_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01


#HIS
en_his <- read.table(file = "Z:/data/mesa_models/MESA.HIS.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
en_his$gene_id <- as.character(en_his$gene_id)
en_his <- subset(en_his, cv_R2_avg > 0.01)
en_his <- en_his[,c(1,10)]
for (i in 1:length(en_his$gene_id)){
  en_his$gene_id[i] <- gsub('\\.[0-9]+','',en_his$gene_id[i])
} #just to remove the decimal places in the gene_id

en_his_2_mets <- read.table(file = "Z:/data/mesa_models/mets_dosages/new_predixcan_HIS2METS_cor", header = T)
en_his_2_mets$gene <- as.character(en_his_2_mets$gene)
filt_en_his_2_mets <- inner_join(en_his, en_his_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01

#AFA
en_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- subset(en_afa, en_afa$cv_R2_avg > 0.01)
en_afa <- en_afa[,c(1,10)]
for (i in 1:length(en_afa$gene_id)){
  en_afa$gene_id[i] <- gsub('\\.[0-9]+','',en_afa$gene_id[i])
} #just to remove the decimal places in the gene_id

en_afa_2_mets <- read.table(file = "Z:/data/mesa_models/mets_dosages/new_predixcan_AFA2METS_cor", header = T)
en_afa_2_mets$gene <- as.character(en_afa_2_mets$gene)
filt_en_afa_2_mets <- inner_join(en_afa, en_afa_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01


library(dplyr)
#RF
#AFA
rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- subset(rf_afa, rf_afa$CV_R2 > 0.01)
for (i in 1:length(rf_afa$Gene_ID)){
  rf_afa$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_afa$Gene_ID[i])
} #just to remove the decimal places in the gene_id

rf_afa_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_rf_cor_test_full_chr.txt", header = T)
rf_afa_2_mets$gene_id <- as.character(rf_afa_2_mets$gene_id)

filt_rf_afa_2_mets <- inner_join(rf_afa, rf_afa_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_rf_afa_2_mets <- filt_rf_afa_2_mets[,c(1,13)]
names(filt_rf_afa_2_mets) <- c("gene", "spearman")

#HIS
rf_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_all_chrom.txt", header = T)
rf_his$Gene_ID <- as.character(rf_his$Gene_ID)
rf_his <- subset(rf_his, rf_his$CV_R2 > 0.01)
for (i in 1:length(rf_his$Gene_ID)){
  rf_his$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_his$Gene_ID[i])
} #just to remove the decimal places in the gene_id

rf_his_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_rf_cor_test_full_chr.txt", header = T)
rf_his_2_mets$gene_id <- as.character(rf_his_2_mets$gene_id)

filt_rf_his_2_mets <- inner_join(rf_his, rf_his_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_rf_his_2_mets <- filt_rf_his_2_mets[,c(1,13)]
names(filt_rf_his_2_mets) <- c("gene", "spearman")

#CAU
rf_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_rf_all_chrom.txt", header = T)
rf_cau$Gene_ID <- as.character(rf_cau$Gene_ID)
rf_cau <- subset(rf_cau, rf_cau$CV_R2 > 0.01)
for (i in 1:length(rf_cau$Gene_ID)){
  rf_cau$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_cau$Gene_ID[i])
} #just to remove the decimal places in the gene_id

rf_cau_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_CAU_2_METS_rf_cor_test_full_chr.txt", header = T)
rf_cau_2_mets$gene_id <- as.character(rf_cau_2_mets$gene_id)

filt_rf_cau_2_mets <- inner_join(rf_cau, rf_cau_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_rf_cau_2_mets <- filt_rf_cau_2_mets[,c(1,13)]
names(filt_rf_cau_2_mets) <- c("gene", "spearman")

#AFHI
rf_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_rf_all_chrom.txt", header = T)
rf_afhi$Gene_ID <- as.character(rf_afhi$Gene_ID)
rf_afhi <- subset(rf_afhi, rf_afhi$CV_R2 > 0.01)
for (i in 1:length(rf_afhi$Gene_ID)){
  rf_afhi$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_afhi$Gene_ID[i])
} #just to remove the decimal places in the gene_id

rf_afhi_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFHI_2_METS_rf_cor_test_full_chr.txt", header = T)
rf_afhi_2_mets$gene_id <- as.character(rf_afhi_2_mets$gene_id)

filt_rf_afhi_2_mets <- inner_join(rf_afhi, rf_afhi_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_rf_afhi_2_mets <- filt_rf_afhi_2_mets[,c(1,13)]
names(filt_rf_afhi_2_mets) <- c("gene", "spearman")


#ALL
rf_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header = T)
rf_all$Gene_ID <- as.character(rf_all$Gene_ID)
rf_all <- subset(rf_all, rf_all$CV_R2 > 0.01)
for (i in 1:length(rf_all$Gene_ID)){
  rf_all$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_all$Gene_ID[i])
} #just to remove the decimal places in the gene_id

rf_all_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_ALL_2_METS_rf_cor_test_full_chr.txt", header = T)
rf_all_2_mets$gene_id <- as.character(rf_all_2_mets$gene_id)

filt_rf_all_2_mets <- inner_join(rf_all, rf_all_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_rf_all_2_mets <- filt_rf_all_2_mets[,c(1,13)]
names(filt_rf_all_2_mets) <- c("gene", "spearman")



#SVR
library(dplyr)
#AFA
svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
svr_afa <- subset(svr_afa, svr_afa$CV_R2 > 0.01)
for (i in 1:length(svr_afa$Gene_ID)){
  svr_afa$Gene_ID[i] <- gsub('\\.[0-9]+','',svr_afa$Gene_ID[i])
} #just to remove the decimal places in the gene_id

svr_afa_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_svr_cor_test_full_chr.txt", header = T)
svr_afa_2_mets$gene_id <- as.character(svr_afa_2_mets$gene_id)

filt_svr_afa_2_mets <- inner_join(svr_afa, svr_afa_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_svr_afa_2_mets <- filt_svr_afa_2_mets[,c(1,13)]
names(filt_svr_afa_2_mets) <- c("gene", "spearman")

#HIS
svr_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_svr_all_chrom.txt", header = T)
svr_his$Gene_ID <- as.character(svr_his$Gene_ID)
svr_his <- subset(svr_his, svr_his$CV_R2 > 0.01)
for (i in 1:length(svr_his$Gene_ID)){
  svr_his$Gene_ID[i] <- gsub('\\.[0-9]+','',svr_his$Gene_ID[i])
} #just to remove the decimal places in the gene_id

svr_his_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_svr_cor_test_full_chr.txt", header = T)
svr_his_2_mets$gene_id <- as.character(svr_his_2_mets$gene_id)

filt_svr_his_2_mets <- inner_join(svr_his, svr_his_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_svr_his_2_mets <- filt_svr_his_2_mets[,c(1,13)]
names(filt_svr_his_2_mets) <- c("gene", "spearman")


#CAU
svr_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_svr_all_chrom.txt", header = T)
svr_cau$Gene_ID <- as.character(svr_cau$Gene_ID)
svr_cau <- subset(svr_cau, svr_cau$CV_R2 > 0.01)
for (i in 1:length(svr_cau$Gene_ID)){
  svr_cau$Gene_ID[i] <- gsub('\\.[0-9]+','',svr_cau$Gene_ID[i])
} #just to remove the decimal places in the gene_id

svr_cau_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_CAU_2_METS_svr_cor_test_full_chr.txt", header = T)
svr_cau_2_mets$gene_id <- as.character(svr_cau_2_mets$gene_id)

filt_svr_cau_2_mets <- inner_join(svr_cau, svr_cau_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_svr_cau_2_mets <- filt_svr_cau_2_mets[,c(1,13)]
names(filt_svr_cau_2_mets) <- c("gene", "spearman")

#AFHI
svr_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_svr_all_chrom.txt", header = T)
svr_afhi$Gene_ID <- as.character(svr_afhi$Gene_ID)
svr_afhi <- subset(svr_afhi, svr_afhi$CV_R2 > 0.01)
for (i in 1:length(svr_afhi$Gene_ID)){
  svr_afhi$Gene_ID[i] <- gsub('\\.[0-9]+','',svr_afhi$Gene_ID[i])
} #just to remove the decimal places in the gene_id

svr_afhi_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFHI_2_METS_svr_cor_test_full_chr.txt", header = T)
svr_afhi_2_mets$gene_id <- as.character(svr_afhi_2_mets$gene_id)

filt_svr_afhi_2_mets <- inner_join(svr_afhi, svr_afhi_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_svr_afhi_2_mets <- filt_svr_afhi_2_mets[,c(1,13)]
names(filt_svr_afhi_2_mets) <- c("gene", "spearman")

#ALL
svr_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header = T)
svr_all$Gene_ID <- as.character(svr_all$Gene_ID)
svr_all <- subset(svr_all, svr_all$CV_R2 > 0.01)
for (i in 1:length(svr_all$Gene_ID)){
  svr_all$Gene_ID[i] <- gsub('\\.[0-9]+','',svr_all$Gene_ID[i])
} #just to remove the decimal places in the gene_id

svr_all_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_ALL_2_METS_svr_cor_test_full_chr.txt", header = T)
svr_all_2_mets$gene_id <- as.character(svr_all_2_mets$gene_id)

filt_svr_all_2_mets <- inner_join(svr_all, svr_all_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_svr_all_2_mets <- filt_svr_all_2_mets[,c(1,13)]
names(filt_svr_all_2_mets) <- c("gene", "spearman")



#KNN
library(dplyr)
#AFA
knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
knn_afa <- subset(knn_afa, knn_afa$CV_R2 > 0.01)
for (i in 1:length(knn_afa$Gene_ID)){
  knn_afa$Gene_ID[i] <- gsub('\\.[0-9]+','',knn_afa$Gene_ID[i])
} #just to remove the decimal places in the gene_id

knn_afa_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_knn_cor_test_full_chr.txt", header = T)
knn_afa_2_mets$gene_id <- as.character(knn_afa_2_mets$gene_id)

filt_knn_afa_2_mets <- inner_join(knn_afa, knn_afa_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_knn_afa_2_mets <- filt_knn_afa_2_mets[,c(1,13)]
names(filt_knn_afa_2_mets) <- c("gene", "spearman")

#HIS
knn_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_knn_all_chrom.txt", header = T)
knn_his$Gene_ID <- as.character(knn_his$Gene_ID)
knn_his <- subset(knn_his, knn_his$CV_R2 > 0.01)
for (i in 1:length(knn_his$Gene_ID)){
  knn_his$Gene_ID[i] <- gsub('\\.[0-9]+','',knn_his$Gene_ID[i])
} #just to remove the decimal places in the gene_id

knn_his_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_knn_cor_test_full_chr.txt", header = T)
knn_his_2_mets$gene_id <- as.character(knn_his_2_mets$gene_id)

filt_knn_his_2_mets <- inner_join(knn_his, knn_his_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_knn_his_2_mets <- filt_knn_his_2_mets[,c(1,13)]
names(filt_knn_his_2_mets) <- c("gene", "spearman")


#CAU
knn_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_knn_all_chrom.txt", header = T)
knn_cau$Gene_ID <- as.character(knn_cau$Gene_ID)
knn_cau <- subset(knn_cau, knn_cau$CV_R2 > 0.01)
for (i in 1:length(knn_cau$Gene_ID)){
  knn_cau$Gene_ID[i] <- gsub('\\.[0-9]+','',knn_cau$Gene_ID[i])
} #just to remove the decimal places in the gene_id

knn_cau_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_CAU_2_METS_knn_cor_test_full_chr.txt", header = T)
knn_cau_2_mets$gene_id <- as.character(knn_cau_2_mets$gene_id)

filt_knn_cau_2_mets <- inner_join(knn_cau, knn_cau_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_knn_cau_2_mets <- filt_knn_cau_2_mets[,c(1,13)]
names(filt_knn_cau_2_mets) <- c("gene", "spearman")


#AFHI
knn_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_knn_all_chrom.txt", header = T)
knn_afhi$Gene_ID <- as.character(knn_afhi$Gene_ID)
knn_afhi <- subset(knn_afhi, knn_afhi$CV_R2 > 0.01)
for (i in 1:length(knn_afhi$Gene_ID)){
  knn_afhi$Gene_ID[i] <- gsub('\\.[0-9]+','',knn_afhi$Gene_ID[i])
} #just to remove the decimal places in the gene_id

knn_afhi_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFHI_2_METS_knn_cor_test_full_chr.txt", header = T)
knn_afhi_2_mets$gene_id <- as.character(knn_afhi_2_mets$gene_id)

filt_knn_afhi_2_mets <- inner_join(knn_afhi, knn_afhi_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_knn_afhi_2_mets <- filt_knn_afhi_2_mets[,c(1,13)]
names(filt_knn_afhi_2_mets) <- c("gene", "spearman")

#ALL
knn_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header = T)
knn_all$Gene_ID <- as.character(knn_all$Gene_ID)
knn_all <- subset(knn_all, knn_all$CV_R2 > 0.01)
for (i in 1:length(knn_all$Gene_ID)){
  knn_all$Gene_ID[i] <- gsub('\\.[0-9]+','',knn_all$Gene_ID[i])
} #just to remove the decimal places in the gene_id

knn_all_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_ALL_2_METS_knn_cor_test_full_chr.txt", header = T)
knn_all_2_mets$gene_id <- as.character(knn_all_2_mets$gene_id)

filt_knn_all_2_mets <- inner_join(knn_all, knn_all_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_knn_all_2_mets <- filt_knn_all_2_mets[,c(1,13)]
names(filt_knn_all_2_mets) <- c("gene", "spearman")


#############make violin plot of all data in table 2, and also do the model comparison for all the pops and AFA and AFHI
#Figure 2
en_afa <- data.frame(spearman=filt_en_afa_2_mets$spearman,Model=rep("EN",nrow(filt_en_afa_2_mets)), 
                     mesa=rep("AFA", nrow(filt_en_afa_2_mets)))
en_his <- data.frame(spearman=filt_en_his_2_mets$spearman,Model=rep("EN",nrow(filt_en_his_2_mets)), 
                     mesa=rep("HIS", nrow(filt_en_his_2_mets)))
en_cau <- data.frame(spearman=filt_en_cau_2_mets$spearman,Model=rep("EN",nrow(filt_en_cau_2_mets)), 
                     mesa=rep("CAU", nrow(filt_en_cau_2_mets)))
en_afhi <- data.frame(spearman=filt_en_afhi_2_mets$spearman,Model=rep("EN",nrow(filt_en_afhi_2_mets)), 
                      mesa=rep("AFHI", nrow(filt_en_afhi_2_mets)))
en_all <- data.frame(spearman=filt_en_all_2_mets$spearman,Model=rep("EN",nrow(filt_en_all_2_mets)), 
                     mesa=rep("ALL", nrow(filt_en_all_2_mets)))


rf_afa <- data.frame(spearman=filt_rf_afa_2_mets$spearman,Model=rep("RF",nrow(filt_rf_afa_2_mets)), 
                     mesa=rep("AFA", nrow(filt_rf_afa_2_mets)))
rf_his <- data.frame(spearman=filt_rf_his_2_mets$spearman,Model=rep("RF",nrow(filt_rf_his_2_mets)), 
                     mesa=rep("HIS", nrow(filt_rf_his_2_mets)))
rf_cau <- data.frame(spearman=filt_rf_cau_2_mets$spearman,Model=rep("RF",nrow(filt_rf_cau_2_mets)), 
                     mesa=rep("CAU", nrow(filt_rf_cau_2_mets)))
rf_afhi <- data.frame(spearman=filt_rf_afhi_2_mets$spearman,Model=rep("RF",nrow(filt_rf_afhi_2_mets)), 
                      mesa=rep("AFHI", nrow(filt_rf_afhi_2_mets)))
rf_all <- data.frame(spearman=filt_rf_all_2_mets$spearman,Model=rep("RF",nrow(filt_rf_all_2_mets)), 
                     mesa=rep("ALL", nrow(filt_rf_all_2_mets)))

svr_afa <- data.frame(spearman=filt_svr_afa_2_mets$spearman,Model=rep("SVR",nrow(filt_svr_afa_2_mets)), 
                      mesa=rep("AFA", nrow(filt_svr_afa_2_mets)))
svr_his <- data.frame(spearman=filt_svr_his_2_mets$spearman,Model=rep("SVR",nrow(filt_svr_his_2_mets)), 
                      mesa=rep("HIS", nrow(filt_svr_his_2_mets)))
svr_cau <- data.frame(spearman=filt_svr_cau_2_mets$spearman,Model=rep("SVR",nrow(filt_svr_cau_2_mets)), 
                      mesa=rep("CAU", nrow(filt_svr_cau_2_mets)))
svr_afhi <- data.frame(spearman=filt_svr_afhi_2_mets$spearman,Model=rep("SVR",nrow(filt_svr_afhi_2_mets)), 
                       mesa=rep("AFHI", nrow(filt_svr_afhi_2_mets)))
svr_all <- data.frame(spearman=filt_svr_all_2_mets$spearman,Model=rep("SVR",nrow(filt_svr_all_2_mets)), 
                      mesa=rep("ALL", nrow(filt_svr_all_2_mets)))

knn_afa <- data.frame(spearman=filt_knn_afa_2_mets$spearman,Model=rep("KNN",nrow(filt_knn_afa_2_mets)), 
                      mesa=rep("AFA", nrow(filt_knn_afa_2_mets)))
knn_his <- data.frame(spearman=filt_knn_his_2_mets$spearman,Model=rep("KNN",nrow(filt_knn_his_2_mets)), 
                      mesa=rep("HIS", nrow(filt_knn_his_2_mets)))
knn_cau <- data.frame(spearman=filt_knn_cau_2_mets$spearman,Model=rep("KNN",nrow(filt_knn_cau_2_mets)), 
                      mesa=rep("CAU", nrow(filt_knn_cau_2_mets)))
knn_afhi <- data.frame(spearman=filt_knn_afhi_2_mets$spearman,Model=rep("KNN",nrow(filt_knn_afhi_2_mets)), 
                       mesa=rep("AFHI", nrow(filt_knn_afhi_2_mets)))
knn_all <- data.frame(spearman=filt_knn_all_2_mets$spearman,Model=rep("KNN",nrow(filt_knn_all_2_mets)), 
                      mesa=rep("ALL", nrow(filt_knn_all_2_mets)))

mesa2mets <- rbind(en_afa, en_his, en_cau, en_afhi, en_all, rf_afa, rf_his, rf_cau, rf_afhi, rf_all,
                   svr_afa, svr_his, svr_cau, svr_afhi, svr_all, knn_afa, knn_his, knn_cau, knn_afhi, knn_all)

mesa2mets_0.1 <- subset(mesa2mets, spearman > 0.1)

library(ggplot2)
ggplot(mesa2mets_0.1, aes(x=mesa, y=spearman, color=Model, fill=Model)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.01, color="black") + theme_classic(20) #+
#scale_y_continuous(breaks=seq(-1.0, 1.0, 0.2), limits=c(-1.0, 1.0)) +
#theme(panel.grid.minor = element_blank()) + xlab("KNN")

ggplot(mesa2mets, aes(x=mesa, y=spearman, color=Model, fill=Model)) + geom_boxplot() + theme_classic(20) +
  xlab("Population") + scale_y_continuous(breaks=seq(-1.0, 1.0, 0.5), limits=c(-1.0, 1.0)) +
  ylab("Spearman Correlation")#stat_summary(fun.y=mean, geom = "point", size=1, color="white")

ggplot(mesa2mets_0.1, aes(x=mesa, y=spearman, color=Model, fill=Model)) + geom_boxplot() + theme_classic(20) +
  xlab("Population") + scale_y_continuous(breaks=seq(0, 1.0, 0.25), limits=c(0, 1.0)) +
  ylab("Spearman Correlation")#stat_summary(fun.y=mean, geom = "point", size=1, color="white")


ggplot(mesa2mets, aes(x=mesa, y=spearman, fill=Model)) + geom_boxplot() + theme_classic(30) +
  xlab("Population") + scale_y_continuous(breaks=seq(-1.0, 1.0, 0.25), limits=c(-1.0, 1.0)) +
  ylab("Spearman Correlation")#stat_summary(fun.y=mean, geom = "point", size=1, color="white")
#width=1000, height=750


#compare ml and en across pops 
#Figure 3

all_en_rf <- inner_join(filt_en_all_2_mets, filt_rf_all_2_mets, by =c("gene_id"="gene"))
all_en_rf <- subset(all_en_rf, spearman.x > 0.1 & spearman.y > 0.1)
all_en_rf <- all_en_rf[,c(3,4)]
names(all_en_rf) <- c("x","y")

p1 <- ggplot(all_en_rf, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70) + labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

all_en_svr <- inner_join(filt_en_all_2_mets, filt_svr_all_2_mets, by =c("gene_id"="gene"))
all_en_svr <- subset(all_en_svr, spearman.x > 0.1 & spearman.y > 0.1)
all_en_svr <- all_en_svr[,c(3,4)]
names(all_en_svr) <- c("x","y")

p2 <- ggplot(all_en_svr, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70) + labs(title="B") + labs(title = "B")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

all_en_knn <- inner_join(filt_en_all_2_mets, filt_knn_all_2_mets, by =c("gene_id"="gene"))
all_en_knn <- subset(all_en_knn, spearman.x > 0.1 & spearman.y > 0.1)
all_en_knn <- all_en_knn[,c(3,4)]
names(all_en_knn) <- c("x","y")

p3 <- ggplot(all_en_knn, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("K Nearest Neighbor") +
  theme_classic(70) + labs(title="C")# theme_classic(20)+ ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

en <- data.frame(spearman=filt_en_all_2_mets$spearman, Model=rep("EN",nrow(filt_en_all_2_mets)))
rf <- data.frame(spearman=filt_rf_all_2_mets$spearman, Model=rep("RF",nrow(filt_rf_all_2_mets)))
svr <- data.frame(spearman=filt_svr_all_2_mets$spearman, Model=rep("SVR",nrow(filt_svr_all_2_mets)))
knn <- data.frame(spearman=filt_knn_all_2_mets$spearman, Model=rep("KNN",nrow(filt_knn_all_2_mets)))

model <- rbind(en, rf, svr, knn)
model_0.1 <- subset(model, spearman > 0.1)

ggplot(model, aes(x = spearman, color=Model)) + 
  scale_x_continuous(name = "Spearman Correlation") +
  geom_density(size=1.5)+ theme_classic(20) + scale_x_continuous(breaks=seq(-1, 1.0, 0.25), limits=c(-1, 1.0)) # +
#geom_vline(data=mu, aes(xintercept=grp.median, color=prediction),linetype="longdash", lwd=1) #+
#scale_color_manual(values = c("red","blue","orange","violet"))

p4 <- ggplot(model_0.1, aes(x = spearman, color=Model)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") +
  geom_density(size=1.5)+ theme_classic(70) + labs(title="D") # +


#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=4000, height=3400

#######
#Make supplementary figure for Figure 3

#Supplementary of Figure 3 
#AFA

afa_en_rf <- inner_join(filt_en_afa_2_mets, filt_rf_afa_2_mets, by =c("gene_id"="gene"))
afa_en_rf <- subset(afa_en_rf, spearman.x > 0.1 & spearman.y > 0.1)
afa_en_rf <- afa_en_rf[,c(3,4)]
names(afa_en_rf) <- c("x","y")

p1 <- ggplot(afa_en_rf, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70) + labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

afa_en_svr <- inner_join(filt_en_afa_2_mets, filt_svr_afa_2_mets, by =c("gene_id"="gene"))
afa_en_svr <- subset(afa_en_svr, spearman.x > 0.1 & spearman.y > 0.1)
afa_en_svr <- afa_en_svr[,c(3,4)]
names(afa_en_svr) <- c("x","y")

p2 <- ggplot(afa_en_svr, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70) + labs(title="B") + labs(title = "B")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

afa_en_knn <- inner_join(filt_en_afa_2_mets, filt_knn_afa_2_mets, by =c("gene_id"="gene"))
afa_en_knn <- subset(afa_en_knn, spearman.x > 0.1 & spearman.y > 0.1)
afa_en_knn <- afa_en_knn[,c(3,4)]
names(afa_en_knn) <- c("x","y")

p3 <- ggplot(afa_en_knn, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("K Nearest Neighbor") +
  theme_classic(70) + labs(title="C")# theme_classic(20)+ ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

en <- data.frame(spearman=filt_en_afa_2_mets$spearman, Model=rep("EN",nrow(filt_en_afa_2_mets)))
rf <- data.frame(spearman=filt_rf_afa_2_mets$spearman, Model=rep("RF",nrow(filt_rf_afa_2_mets)))
svr <- data.frame(spearman=filt_svr_afa_2_mets$spearman, Model=rep("SVR",nrow(filt_svr_afa_2_mets)))
knn <- data.frame(spearman=filt_knn_afa_2_mets$spearman, Model=rep("KNN",nrow(filt_knn_afa_2_mets)))

model <- rbind(en, rf, svr, knn)
model_0.1 <- subset(model, spearman > 0.1)

p4 <- ggplot(model_0.1, aes(x = spearman, color=Model)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") +
  geom_density(size=1.5)+ theme_classic(70) + labs(title="D") # +


#plot afa 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=4000, height=3400


#Supplementary of Figure 3 
#his

his_en_rf <- inner_join(filt_en_his_2_mets, filt_rf_his_2_mets, by =c("gene_id"="gene"))
his_en_rf <- subset(his_en_rf, spearman.x > 0.1 & spearman.y > 0.1)
his_en_rf <- his_en_rf[,c(3,4)]
names(his_en_rf) <- c("x","y")

p1 <- ggplot(his_en_rf, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70) + labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

his_en_svr <- inner_join(filt_en_his_2_mets, filt_svr_his_2_mets, by =c("gene_id"="gene"))
his_en_svr <- subset(his_en_svr, spearman.x > 0.1 & spearman.y > 0.1)
his_en_svr <- his_en_svr[,c(3,4)]
names(his_en_svr) <- c("x","y")

p2 <- ggplot(his_en_svr, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70) + labs(title="B") + labs(title = "B")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

his_en_knn <- inner_join(filt_en_his_2_mets, filt_knn_his_2_mets, by =c("gene_id"="gene"))
his_en_knn <- subset(his_en_knn, spearman.x > 0.1 & spearman.y > 0.1)
his_en_knn <- his_en_knn[,c(3,4)]
names(his_en_knn) <- c("x","y")

p3 <- ggplot(his_en_knn, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("K Nearest Neighbor") +
  theme_classic(70) + labs(title="C")# theme_classic(20)+ ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

en <- data.frame(spearman=filt_en_his_2_mets$spearman, Model=rep("EN",nrow(filt_en_his_2_mets)))
rf <- data.frame(spearman=filt_rf_his_2_mets$spearman, Model=rep("RF",nrow(filt_rf_his_2_mets)))
svr <- data.frame(spearman=filt_svr_his_2_mets$spearman, Model=rep("SVR",nrow(filt_svr_his_2_mets)))
knn <- data.frame(spearman=filt_knn_his_2_mets$spearman, Model=rep("KNN",nrow(filt_knn_his_2_mets)))

model <- rbind(en, rf, svr, knn)
model_0.1 <- subset(model, spearman > 0.1)

p4 <- ggplot(model_0.1, aes(x = spearman, color=Model)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") +
  geom_density(size=1.5)+ theme_classic(70) + labs(title="D") # +


#plot his 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=4000, height=3400



#Supplementary of Figure 3 
#cau

cau_en_rf <- inner_join(filt_en_cau_2_mets, filt_rf_cau_2_mets, by =c("gene_id"="gene"))
cau_en_rf <- subset(cau_en_rf, spearman.x > 0.1 & spearman.y > 0.1)
cau_en_rf <- cau_en_rf[,c(3,4)]
names(cau_en_rf) <- c("x","y")

p1 <- ggplot(cau_en_rf, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70) + labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

cau_en_svr <- inner_join(filt_en_cau_2_mets, filt_svr_cau_2_mets, by =c("gene_id"="gene"))
cau_en_svr <- subset(cau_en_svr, spearman.x > 0.1 & spearman.y > 0.1)
cau_en_svr <- cau_en_svr[,c(3,4)]
names(cau_en_svr) <- c("x","y")

p2 <- ggplot(cau_en_svr, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70) + labs(title="B") + labs(title = "B")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

cau_en_knn <- inner_join(filt_en_cau_2_mets, filt_knn_cau_2_mets, by =c("gene_id"="gene"))
cau_en_knn <- subset(cau_en_knn, spearman.x > 0.1 & spearman.y > 0.1)
cau_en_knn <- cau_en_knn[,c(3,4)]
names(cau_en_knn) <- c("x","y")

p3 <- ggplot(cau_en_knn, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("K Nearest Neighbor") +
  theme_classic(70) + labs(title="C")# theme_classic(20)+ ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

en <- data.frame(spearman=filt_en_cau_2_mets$spearman, Model=rep("EN",nrow(filt_en_cau_2_mets)))
rf <- data.frame(spearman=filt_rf_cau_2_mets$spearman, Model=rep("RF",nrow(filt_rf_cau_2_mets)))
svr <- data.frame(spearman=filt_svr_cau_2_mets$spearman, Model=rep("SVR",nrow(filt_svr_cau_2_mets)))
knn <- data.frame(spearman=filt_knn_cau_2_mets$spearman, Model=rep("KNN",nrow(filt_knn_cau_2_mets)))

model <- rbind(en, rf, svr, knn)
model_0.1 <- subset(model, spearman > 0.1)

p4 <- ggplot(model_0.1, aes(x = spearman, color=Model)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") +
  geom_density(size=1.5)+ theme_classic(70) + labs(title="D") # +


#plot cau 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=4000, height=3400



#Supplementary of Figure 3 
#afhi

afhi_en_rf <- inner_join(filt_en_afhi_2_mets, filt_rf_afhi_2_mets, by =c("gene_id"="gene"))
afhi_en_rf <- subset(afhi_en_rf, spearman.x > 0.1 & spearman.y > 0.1)
afhi_en_rf <- afhi_en_rf[,c(3,4)]
names(afhi_en_rf) <- c("x","y")

p1 <- ggplot(afhi_en_rf, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70) + labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

afhi_en_svr <- inner_join(filt_en_afhi_2_mets, filt_svr_afhi_2_mets, by =c("gene_id"="gene"))
afhi_en_svr <- subset(afhi_en_svr, spearman.x > 0.1 & spearman.y > 0.1)
afhi_en_svr <- afhi_en_svr[,c(3,4)]
names(afhi_en_svr) <- c("x","y")

p2 <- ggplot(afhi_en_svr, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70) + labs(title="B") + labs(title = "B")# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

afhi_en_knn <- inner_join(filt_en_afhi_2_mets, filt_knn_afhi_2_mets, by =c("gene_id"="gene"))
afhi_en_knn <- subset(afhi_en_knn, spearman.x > 0.1 & spearman.y > 0.1)
afhi_en_knn <- afhi_en_knn[,c(3,4)]
names(afhi_en_knn) <- c("x","y")

p3 <- ggplot(afhi_en_knn, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("K Nearest Neighbor") +
  theme_classic(70) + labs(title="C")# theme_classic(20)+ ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

en <- data.frame(spearman=filt_en_afhi_2_mets$spearman, Model=rep("EN",nrow(filt_en_afhi_2_mets)))
rf <- data.frame(spearman=filt_rf_afhi_2_mets$spearman, Model=rep("RF",nrow(filt_rf_afhi_2_mets)))
svr <- data.frame(spearman=filt_svr_afhi_2_mets$spearman, Model=rep("SVR",nrow(filt_svr_afhi_2_mets)))
knn <- data.frame(spearman=filt_knn_afhi_2_mets$spearman, Model=rep("KNN",nrow(filt_knn_afhi_2_mets)))

model <- rbind(en, rf, svr, knn)
model_0.1 <- subset(model, spearman > 0.1)

p4 <- ggplot(model_0.1, aes(x = spearman, color=Model)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") +
  geom_density(size=1.5)+ theme_classic(70) + labs(title="D") # +


#plot afhi 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=4000, height=3400
