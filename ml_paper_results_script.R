# Our ML Performance comparison to Elastic Net Paper
#Results Sections

#Cross Validated Performance of all ML models and Population
#Count genes for R2 > 0.01, 0.05, 0.1, 0.5

#AFA
en_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)

mean(en_afa$cv_R2_avg)
median(en_afa$cv_R2_avg)
en_afa <- subset(en_afa, en_afa$cv_R2_avg > -1)
en_afa <- en_afa[,c(1,10)]
names(en_afa) <- c("gene", "en")

rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
mean(rf_afa$CV_R2)
median(rf_afa$CV_R2)
rf_afa <- subset(rf_afa, rf_afa$CV_R2 > -1)
rf_afa <- rf_afa[,c(1,3)]
names(rf_afa) <- c("gene", "rf")

svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
mean(svr_afa$CV_R2)
median(svr_afa$CV_R2)
svr_afa <- subset(svr_afa, svr_afa$CV_R2 > -1)
svr_afa <- svr_afa[,c(1,3)]
names(svr_afa) <- c("gene", "svr")

knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
mean(knn_afa$CV_R2)
median(knn_afa$CV_R2)
knn_afa <- subset(knn_afa, knn_afa$CV_R2 > -1)
knn_afa <- knn_afa[,c(1,3)]
names(knn_afa) <- c("gene", "knn")


# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_afa, rf_afa, by = c("gene"="gene"))
cor.test(en_rf$en, en_rf$rf)
t.test(en_rf$en, en_rf$rf)
library("ggpubr")
ggscatter(en_rf, x = "en", y = "rf", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "Random Forest", 
          title = "Cross Validated Performance",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_afa, svr_afa, by = c("gene"="gene"))
cor.test(en_svr$en, en_svr$svr)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_afa, knn_afa, by = c("gene"="gene"))
cor.test(en_knn$en, en_knn$knn)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)


#HIS
en_his <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_HIS_model_summaries.txt", header = TRUE)
en_his$gene_id <- as.character(en_his$gene_id)
library(tidyverse)
#en_his <- en_his[,c(2,10)]
en_his <- drop_na(en_his)
mean(en_his$cv_R2_avg)
median(en_his$cv_R2_avg)
#check for NA's in the avg_cv_R2 and drop it
en_his <- subset(en_his, en_his$cv_R2_avg > -1)
en_his <- en_his[,c(1,10)]
names(en_his) <- c("gene","en")

rf_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_all_chrom.txt", header = T)
rf_his$Gene_ID <- as.character(rf_his$Gene_ID)
mean(rf_his$CV_R2)
median(rf_his$CV_R2)
rf_his <- subset(rf_his, rf_his$CV_R2 > -1)
rf_his <- rf_his[,c(1,3)]
names(rf_his) <- c("gene", "rf")

svr_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_svr_all_chrom.txt", header = T)
svr_his$Gene_ID <- as.character(svr_his$Gene_ID)
mean(svr_his$CV_R2)
median(svr_his$CV_R2)
svr_his <- subset(svr_his, svr_his$CV_R2 > -1)
svr_his <- svr_his[,c(1,3)]
names(svr_his) <- c("gene","svr")

knn_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_knn_all_chrom.txt", header = T)
knn_his$Gene_ID <- as.character(knn_his$Gene_ID)
mean(knn_his$CV_R2)
median(knn_his$CV_R2)
knn_his <- subset(knn_his, knn_his$CV_R2 > -1)
knn_his <- knn_his[,c(1,3)]
names(knn_his) <- c("gene", "knn")

# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_his, rf_his, by = c("gene"="gene"))
cor.test(en_rf$en, en_rf$rf)
t.test(en_rf$en, en_rf$rf)

library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_his, svr_his, by = c("gene"="gene"))
cor.test(en_svr$en, en_svr$svr)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_his, knn_his, by = c("gene"="gene"))
cor.test(en_knn$en, en_knn$knn)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)




#CAU
en_cau <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_CAU_model_summaries.txt", header = TRUE)
en_cau$gene_id <- as.character(en_cau$gene_id)
library(tidyverse)
en_cau <- drop_na(en_cau)
mean(en_cau$cv_R2_avg)
median(en_cau$cv_R2_avg)
#check for NA's in the avg_cv_R2 and drop it
en_cau <- subset(en_cau, en_cau$cv_R2_avg > -1)
en_cau <- en_cau[,c(1,10)]
names(en_cau) <- c("gene","en")

rf_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_rf_all_chrom.txt", header = T)
rf_cau$Gene_ID <- as.character(rf_cau$Gene_ID)
mean(rf_cau$CV_R2)
median(rf_cau$CV_R2)
rf_cau <- subset(rf_cau, rf_cau$CV_R2 > -1)
rf_cau <- rf_cau[,c(1,3)]
names(rf_cau) <- c("gene", "rf")

svr_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_svr_all_chrom.txt", header = T)
svr_cau$Gene_ID <- as.character(svr_cau$Gene_ID)
mean(svr_cau$CV_R2)
median(svr_cau$CV_R2)
svr_cau <- subset(svr_cau, svr_cau$CV_R2 > -1)
svr_cau <- svr_cau[,c(1,3)]
names(svr_cau) <- c("gene","svr")

knn_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_knn_all_chrom.txt", header = T)
knn_cau$Gene_ID <- as.character(knn_cau$Gene_ID)
mean(knn_cau$CV_R2)
median(knn_cau$CV_R2)
knn_cau <- subset(knn_cau, knn_cau$CV_R2 > -1)
knn_cau <- knn_cau[,c(1,3)]
names(knn_cau) <- c("gene", "knn")

# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_cau, rf_cau, by = c("gene"="gene"))
cor.test(en_rf$en, en_rf$rf)
t.test(en_rf$en, en_rf$rf)

library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_cau, svr_cau, by = c("gene"="gene"))
cor.test(en_svr$en, en_svr$svr)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_cau, knn_cau, by = c("gene"="gene"))
cor.test(en_knn$en, en_knn$knn)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)


####################################################
#check if the genes in the models after filtering are the same
library(dplyr)
knng <- semi_join(knn_his, en_his, by=c("Gene_ID" = "gene_id"))
knng2 <- anti_join(knn_his, en_his, by=c("Gene_ID" = "gene_id"))

svrg <- semi_join(svr_his, en_his, by=c("Gene_ID" = "gene_id"))
svrg2 <- anti_join(svr_his, en_his, by=c("Gene_ID" = "gene_id"))

rfg <- semi_join(rf_his, en_his, by=c("Gene_ID" = "gene_id"))
rfg2 <- anti_join(rf_his, en_his, by=c("Gene_ID" = "gene_id"))

eng <- semi_join(en_his, rf_his, by=c("gene_id" = "Gene_ID"))
eng2 <- anti_join(en_his, rf_his, by=c("gene_id" = "Gene_ID"))

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




#compare the R2 of HLA-DRB1 and HLA-DRB5 for the four models
r2 <- c(0.7985787,0.8731306,0.7552526,0.6289000, 0.8098494, 0.8815585, 0.7609722, 0.6093821)
d1 <- c(0.7985787,0.8731306,0.7552526,0.6289000)
d1 <- data.frame(d1)
d1$group <- c("EN", "RF", "SVR", "KNN")
res.aov <- aov(d1 ~ group, data=d1)
summary(res.aov)

drb_df <- matrix(data=r2, ncol=2)
row.names(drb_df) <- c("EN", "RF", "SVR", "KNN")
colnames(drb_df) <- c("HLA-DRB1", "HLA-DRB5")



#