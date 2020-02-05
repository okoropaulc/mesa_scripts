library(dplyr)
unkn <- anti_join(rf_all, en_all, by = c("Gene_ID"="gene_id"))

en0.8 <- subset(en_all, cv_R2_avg > 0.8)
rf0.8 <- subset(rf_all, CV_R2 > 0.8)
svr0.8 <- subset(svr_all, CV_R2 > 0.8)

j08 <- inner_join(en0.8, rf0.8, by=c("gene_id"="Gene_ID"))
j08 <- inner_join(j08, svr0.8, by=c("gene_id"="Gene_ID"))

enf <- subset(en_all, cv_R2_avg >= 0.5)
rff <- subset(rf_all, CV_R2 >= 0.5)
svrf <- subset(svr_all, CV_R2 >= 0.5)
knnf <- subset(knn_all, CV_R2 >= 0.5)

jf <- inner_join(enf, rff, by=c("gene_id"="Gene_ID"))
jf <- inner_join(jf, svrf, by=c("gene_id"="Gene_ID"))

jf$gene_name <- as.character(jf$gene_name)

#do fuma
rf_all$Gene_Name <- as.character(rf_all$Gene_Name)
write.table(rf_all$Gene_Name, file="C:/Users/okoro/OneDrive/Desktop/all_genes.txt",row.names=F,col.names=F,quote=F)
write.table(jf$gene_name, file="C:/Users/okoro/OneDrive/Desktop/all0.5_genes.txt",row.names=F,col.names=F,quote=F)


#filter the mesa models with R2 and filter the spearman
library(dplyr)
#CAU
en_cau <- read.table(file = "Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
en_cau$gene_id <- as.character(en_cau$gene_id)
en_cau <- subset(en_cau, en_cau$cv_R2_avg > 0.01)
en_cau <- en_cau[,c(1,10)]
for (i in 1:length(en_cau$gene_id)){
  en_cau$gene_id[i] <- gsub('\\.[0-9]+','',en_cau$gene_id[i])
} #just to remove the decimal places in the gene_id

en_cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
en_cau_2_mets$gene <- as.character(en_cau_2_mets$gene)
filt_en_cau_2_mets <- inner_join(en_cau, en_cau_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01
mean(filt_en_cau_2_mets$spearman)

#HIS
en_his <- read.table(file = "Z:/data/mesa_models/MESA.HIS.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
en_his$gene_id <- as.character(en_his$gene_id)
en_his <- subset(en_his, cv_R2_avg > 0.01)
en_his <- en_his[,c(1,10)]
for (i in 1:length(en_his$gene_id)){
  en_his$gene_id[i] <- gsub('\\.[0-9]+','',en_his$gene_id[i])
} #just to remove the decimal places in the gene_id

en_his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
en_his_2_mets$gene <- as.character(en_his_2_mets$gene)
filt_en_his_2_mets <- inner_join(en_his, en_his_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01
mean(filt_en_his_2_mets$spearman)


#AFA
en_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- subset(en_afa, en_afa$cv_R2_avg > 0.01)
en_afa <- en_afa[,c(1,10)]
for (i in 1:length(en_afa$gene_id)){
  en_afa$gene_id[i] <- gsub('\\.[0-9]+','',en_afa$gene_id[i])
} #just to remove the decimal places in the gene_id

en_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
en_afa_2_mets$gene <- as.character(en_afa_2_mets$gene)
filt_en_afa_2_mets <- inner_join(en_afa, en_afa_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01
mean(filt_en_afa_2_mets$spearman)


