# Our ML Performance comparison to Elastic Net Paper
#Results Sections

#Cross Validated Performance of all ML models and Population
#Count genes for R2 > 0.01, 0.05, 0.1, 0.5

#AFA
en_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- subset(en_afa, en_afa$cv_R2_avg > 0.5)

rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- subset(rf_afa, rf_afa$CV_R2 > 0.5)

svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
svr_afa <- subset(svr_afa, svr_afa$CV_R2 > 0.5)

knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
knn_afa <- subset(knn_afa, knn_afa$CV_R2 > 0.01)


#HIS
en_his <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_HIS_model_summaries.txt", header = TRUE)
en_his$gene_id <- as.character(en_his$gene_id)
#check for NA's in the avg_cv_R2 and drop it
library(tidyverse)
en_his <- en_his[,c(2,10)]
en_his <- drop_na(en_his)
en_his <- subset(en_his, en_his$cv_R2_avg > 0.01)

rf_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_all_chrom.txt", header = T)
#rf_his$Gene_ID <- as.character(rf_his$Gene_ID)
rf_his <- subset(rf_his, rf_his$CV_R2 > 0.01)

svr_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_svr_all_chrom.txt", header = T)
#svr_his$Gene_ID <- as.character(svr_his$Gene_ID)
svr_his <- subset(svr_his, svr_his$CV_R2 > 0.01)

knn_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_knn_all_chrom.txt", header = T)
knn_his$Gene_ID <- as.character(knn_his$Gene_ID)
knn_his <- subset(knn_his, knn_his$CV_R2 > 0.01)



library(dplyr)


############
#Do FUMA for Genes in EN and RF with R > 0.3 AFA 2 METS
elnet_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
elnet_0.3 <- subset(elnet_afa_2_mets, spearman > 0.3)

#RF
rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- subset(rf_afa, rf_afa$CV_R2 > 0.01)
for (i in 1:length(rf_afa$Gene_ID)){
  rf_afa$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_afa$Gene_ID[i])
} #just to remove the decimal places in the gene_id


rf_afa_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_rf_cor_test_full_chr.txt", header = T)
rf_afa_2_mets$gene_id <- as.character(rf_afa_2_mets$gene_id)

filt_rf_afa_2_mets <- inner_join(rf_afa, rf_afa_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_rf_afa_2_mets <- filt_rf_afa_2_mets[,c(1,13)]
names(filt_rf_afa_2_mets) <- c("gene", "rf_spearman")

afa_rf_0.3 <- subset(filt_rf_afa_2_mets, rf_spearman > 0.3)

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
gv18$gene_id <- as.character(gv18$gene_id)
gv18_intact <- gv18
for (i in 1:length(gv18$gene_id)){
  gv18$gene_id[i] <- gsub('\\.[0-9]+','',gv18$gene_id[i])
} #just to remove the decimal places in the gene_id

gv18 <- subset(gv18, gene_type == "protein_coding")

#All 9623 genes used for model building
rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- subset(rf_afa, rf_afa$CV_R2 > 0.01)
for (i in 1:length(rf_afa$Gene_ID)){
  rf_afa$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_afa$Gene_ID[i])
} #just to remove the decimal places in the gene_id

write.table(rf_afa[2], file = "/Users/okoro/OneDrive/Desktop/elnet_all.txt", quote = FALSE, row.names = FALSE)

write.table(afa_rf_0.3[1], file = "/Users/okoro/OneDrive/Desktop/rf_afa_0.3.txt", quote = FALSE, row.names = FALSE)


#for genes in EN or RF or Intersect R > 0.1
elnet_afa_2_mets_0.1 <- subset(elnet_afa_2_mets, spearman > 0.1)
filt_rf_afa_2_mets_0.1 <- subset(filt_rf_afa_2_mets, rf_spearman > 0.1)
en_rf_afa_intersect <- inner_join(elnet_afa_2_mets_0.1, filt_rf_afa_2_mets_0.1, by = c("gene" = "gene"))
write.table(en_rf_afa_intersect[1], file = "/Users/okoro/OneDrive/Desktop/en_rf_afa_inter.txt", quote = FALSE, row.names = FALSE)

filt_elnet_rf_afa_2_mets <- inner_join(elnet_afa_2_mets, filt_rf_afa_2_mets, by = c("gene" = "gene"))
sub_elnet_rf_afa_2_mets <- subset(filt_elnet_rf_afa_2_mets, spearman > 0.1 | rf_spearman > 0.1) #removes where both are < 0.1

rfonly <- anti_join(filt_rf_afa_2_mets_0.1, elnet_afa_2_mets_0.1, by = c("gene" = "gene")) #Better!
write.table(rfonly[1], file = "/Users/okoro/OneDrive/Desktop/rf_only.txt", quote = FALSE, row.names = FALSE)

enonly <- anti_join(elnet_afa_2_mets_0.1, filt_rf_afa_2_mets_0.1,by = c("gene" = "gene"))
write.table(enonly[1], file = "/Users/okoro/OneDrive/Desktop/en_only.txt", quote = FALSE, row.names = FALSE)
