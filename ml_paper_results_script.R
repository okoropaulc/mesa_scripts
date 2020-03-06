# Our ML Performance comparison to Elastic Net Paper
#Results Sections

#Cross Validated Performance of all ML models and Population
#Count genes for R2 > 0.01, 0.05, 0.1, 0.5
#
#the new_mesa folder in ryan directory is what i used for cv R2. Even though the model summaries are named different,
#they are same files

#CAU
en_cau <- read.table(file = "Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
en_cau <- subset(en_cau, cv_R2_avg > 0.5)
en_cau <- en_cau[,c(1,10)]
names(en_cau) <- c("gene", "en")

#en_cau2 <- read.table(file = "Z:/data/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
#en_cau2 <- subset(en_cau2, cv_R2_avg > 0.5)
#HIS
en_his <- read.table(file = "Z:/data/mesa_models/MESA.HIS.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
en_his <- subset(en_his, cv_R2_avg > 0.5)
en_his <- en_his[,c(1,10)]
names(en_his) <- c("gene", "en")

#AFA
en_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)

mean(en_afa$cv_R2_avg)
median(en_afa$cv_R2_avg)
en_afa <- subset(en_afa, en_afa$cv_R2_avg > 0.5)
en_afa <- en_afa[,c(1,10)]
names(en_afa) <- c("gene", "en")

#best_grid_rf_all_chrom.txt
rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
mean(rf_afa$CV_R2)
median(rf_afa$CV_R2)
rf_afa <- subset(rf_afa, rf_afa$CV_R2 > 0.5)
rf_afa <- rf_afa[,c(1,3)]
names(rf_afa) <- c("gene", "rf")

svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
mean(svr_afa$CV_R2)
median(svr_afa$CV_R2)
svr_afa <- subset(svr_afa, svr_afa$CV_R2 > 0.5)
svr_afa <- svr_afa[,c(1,3)]
names(svr_afa) <- c("gene", "svr")

knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
mean(knn_afa$CV_R2)
median(knn_afa$CV_R2)
knn_afa <- subset(knn_afa, knn_afa$CV_R2 > 0.5)
knn_afa <- knn_afa[,c(1,3)]
names(knn_afa) <- c("gene", "knn")

#check genes across models and populations
library(dplyr)
afa_genes <- inner_join(en_afa, knn_afa, by=c("gene_id"="Gene_ID"))
afa_genes <- inner_join(afa_genes, svr_afa, by=c("gene_id"="Gene_ID"))
afa_genes <- inner_join(afa_genes, rf_afa, by=c("gene_id"="Gene_ID"))

his_genes <- inner_join(en_his, knn_afa, by=c("gene_id"="Gene_ID"))
his_genes <- inner_join(his_genes, svr_afa, by=c("gene_id"="Gene_ID"))
his_genes <- inner_join(his_genes, rf_afa, by=c("gene_id"="Gene_ID"))

cau_genes <- inner_join(en_cau, knn_afa, by=c("gene_id"="Gene_ID"))
cau_genes <- inner_join(cau_genes, svr_afa, by=c("gene_id"="Gene_ID"))
cau_genes <- inner_join(cau_genes, rf_afa, by=c("gene_id"="Gene_ID"))

all_genes <- inner_join(en_all, knn_afa, by=c("gene_id"="Gene_ID"))
all_genes <- inner_join(all_genes, svr_afa, by=c("gene_id"="Gene_ID"))
all_genes <- inner_join(all_genes, rf_afa, by=c("gene_id"="Gene_ID"))

afa_his_cau_genes <- inner_join(afa_genes, his_genes, by=c("gene_id"="gene_id"))
afa_his_cau_genes <- inner_join(afa_his_cau_genes, cau_genes, by=c("gene_id"="gene_id"))

afa_his_cau_all_genes <- inner_join(afa_his_cau_genes, all_genes, by=c("gene_id"="gene_id"))

# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_cau, rf_afa, by = c("gene"="gene"))
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
p1 <- ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(50)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_cau, svr_afa, by = c("gene"="gene"))
cor.test(en_svr$en, en_svr$svr)
library(ggplot2) #i used geom smooth to add the red regression line
p2 <- ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(50)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_cau, knn_afa, by = c("gene"="gene"))
cor.test(en_knn$en, en_knn$knn)
library(ggplot2) #i used geom smooth to add the red regression line
p3 <- ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(50)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)


####################
#plot violin plot of the cv r2 of en, rf, svr, knn
en <- data.frame(alg=rep("EN", nrow(en_cau)), r2=en_cau$en)
rf <- data.frame(alg=rep("RF", nrow(rf_afa)), r2=rf_afa$rf)
svr <- data.frame(alg=rep("SVR", nrow(svr_afa)), r2=svr_afa$svr)
knn <- data.frame(alg=rep("KNN", nrow(knn_afa)), r2=knn_afa$knn)

r2_comp <- rbind(en,rf,svr,knn)
p4 <- ggplot(r2_comp, aes(x=alg, y=r2, color=alg, fill=alg)) + 
  geom_violin(trim=F) + theme_classic(50) + xlab("") + theme(legend.position = "none") +
  geom_boxplot(width=0.15,color="black") + ylab("CV Performance") + xlab("Model")

#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=3500, height=2500




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
en_cau <- subset(en_cau, en_cau$cv_R2_avg > 0.01)
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



####AFHI and ALL summaries need header
en_cau <- read.table(file="Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
#en_all <- fread(file="Z:/no_header_MESA.ALL.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=F)
header <- colnames(en_cau)
#names(en_all) <- header

en_all$gene_id <- as.character(en_all$gene_id)
for (i in 1:nrow(en_all)){
  en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
} #just to remove the decimal places in the gene_id


en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)

#filter by cv R2 > 0.01
en_all <- subset(en_all, cv_R2_avg > 0.01)
#####AFHI
en_afhi <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFHI_model_summaries.txt", header=T)#it has no header
#names(en_afhi) <- header
#write.table(en_afhi, file="Z:/data/mesa_models/split_mesa/results/all_chr_AFHI_model_summaries.txt", row.names=F, quote=F, sep="\t")
en_afhi <- en_afhi[,c(1,2,10)]
en_afhi$cv_R2_avg <- as.numeric(en_afhi$cv_R2_avg)
mean(en_afhi$cv_R2_avg)
median(en_afhi$cv_R2_avg)
en_afhi <- subset(en_afhi, cv_R2_avg > 0.1)
en_afhi$gene_id <- as.character(en_afhi$gene_id)


rf_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_rf_all_chrom.txt", header=T)
rf_afhi$CV_R2 <- as.numeric(rf_afhi$CV_R2)
mean(rf_afhi$CV_R2)
median(rf_afhi$CV_R2)
rf_afhi <- subset(rf_afhi, CV_R2 > 0.1)
rf_afhi$Gene_ID <- as.character(rf_afhi$Gene_ID)


svr_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_svr_all_chrom.txt", header=T)
svr_afhi$CV_R2 <- as.numeric(svr_afhi$CV_R2)
mean(svr_afhi$CV_R2)
median(svr_afhi$CV_R2)
svr_afhi <- subset(svr_afhi, CV_R2 > 0.1)
svr_afhi$Gene_ID <- as.character(svr_afhi$Gene_ID)


knn_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_knn_all_chrom.txt", header=T)
knn_afhi$CV_R2 <- as.numeric(knn_afhi$CV_R2)
mean(knn_afhi$CV_R2)
median(knn_afhi$CV_R2)
knn_afhi <- subset(knn_afhi, CV_R2 > 0.1)
knn_afhi$Gene_ID <- as.character(knn_afhi$Gene_ID)



# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_afhi, rf_afhi, by = c("gene_id"="Gene_ID"))
en_rf <- en_rf[,c(3,5)]
names(en_rf) <- c("en", "rf")
cor.test(en_rf$en, en_rf$rf)
t.test(en_rf$en, en_rf$rf)

library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_afhi, svr_afhi, by = c("gene_id"="Gene_ID"))
en_svr <- en_svr[,c(3,5)]
names(en_svr) <- c("en", "svr")
cor.test(en_svr$en, en_svr$svr)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_afhi, knn_afhi, by = c("gene_id"="Gene_ID"))
en_knn <- en_knn[,c(3,5)]
names(en_knn) <- c("en","knn")
cor.test(en_knn$en, en_knn$knn)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)


######################
###########               ALL

en_all <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", header=T)#it has no header
#names(en_all) <- header
#write.table(en_all,file="Z:/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", row.names=F, quote=F, sep="\t")
en_all <- en_all[,c(1,2,10)]
library(tidyverse)
en_all <- drop_na(en_all) #remove NA
en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)
mean(en_all$cv_R2_avg)
median(en_all$cv_R2_avg)
en_all <- subset(en_all, cv_R2_avg > 0.5)
en_all$gene_id <- as.character(en_all$gene_id)


rf_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
rf_all$CV_R2 <- as.numeric(rf_all$CV_R2)
mean(rf_all$CV_R2)
median(rf_all$CV_R2)
rf_all <- subset(rf_all, CV_R2 > 0.1)
rf_all$Gene_ID <- as.character(rf_all$Gene_ID)


svr_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
svr_all$CV_R2 <- as.numeric(svr_all$CV_R2)
mean(svr_all$CV_R2)
median(svr_all$CV_R2)
svr_all <- subset(svr_all, CV_R2 > 0.1)
svr_all$Gene_ID <- as.character(svr_all$Gene_ID)


knn_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)
knn_all$CV_R2 <- as.numeric(knn_all$CV_R2)
mean(knn_all$CV_R2)
median(knn_all$CV_R2)
knn_all <- subset(knn_all, CV_R2 > 0.1)
knn_all$Gene_ID <- as.character(knn_all$Gene_ID)



# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_all, rf_all, by = c("gene_id"="Gene_ID"))
en_rf <- en_rf[,c(3,5)]
names(en_rf) <- c("en", "rf")
cor.test(en_rf$en, en_rf$rf)
t.test(en_rf$en, en_rf$rf)

library(ggplot2) #i used geom smooth to add the red regression line
p1 <- ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("Random Forest") + theme_classic(30) +
  xlim(-0.1,1) + ylim(-1,1) # + ggtitle("Cross Validation Perfromance") (save at width=900, height=600) + xlim(-1,1) + ylim(-1,1)

en_svr <- inner_join(en_all, svr_all, by = c("gene_id"="Gene_ID"))
en_svr <- en_svr[,c(3,5)]
names(en_svr) <- c("en", "svr")
cor.test(en_svr$en, en_svr$svr)
t.test(en_svr$en, en_svr$svr)
library(ggplot2) #i used geom smooth to add the red regression line
p2 <- ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("Support Vector") + theme_classic(30) + 
  xlim(-0.1,1) + ylim(-1,1) # + ggtitle("Cross Validation Perfromance") (save at width=900, height=600)+ xlim(-1,1) + ylim(-1,1)

en_knn <- inner_join(en_all, knn_all, by = c("gene_id"="Gene_ID"))
en_knn <- en_knn[,c(3,5)]
names(en_knn) <- c("en","knn")
cor.test(en_knn$en, en_knn$knn)
t.test(en_knn$en, en_knn$knn)
library(ggplot2) #i used geom smooth to add the red regression line
p3 <- ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("K-Nearest Neighbor") + theme_classic(30) +
  xlim(-0.1,1) + ylim(-1,1) # + ggtitle("Cross Validation Perfromance") (save at width=850, height=600) + xlim(-1,1) + ylim(-1,1)


####################
#plot violin plot of the cv r2 of en, rf, svr, knn
en <- data.frame(alg=rep("EN", nrow(en_all)), r2=en_all$cv_R2_avg)
rf <- data.frame(alg=rep("RF", nrow(rf_all)), r2=rf_all$CV_R2)
svr <- data.frame(alg=rep("SVR", nrow(svr_all)), r2=svr_all$CV_R2)
knn <- data.frame(alg=rep("KNN", nrow(knn_all)), r2=knn_all$CV_R2)

r2_comp <- rbind(en,rf,svr,knn)
p4 <- ggplot(r2_comp, aes(x=alg, y=r2, color=alg, fill=alg)) + 
  geom_violin(trim=F) + theme_classic(30) + xlab("") + theme(legend.position = "none") +
  geom_boxplot(width=0.15,color="black") + ylab("CV Performance")

#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=3500, height=2500


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
mean(filt_en_afhi_2_mets$spearman)
filt_en_afhi_2_mets <- subset(filt_en_afhi_2_mets, spearman > 0.1)

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
mean(filt_en_all_2_mets$spearman)
filt_en_all_2_mets <- subset(filt_en_all_2_mets, spearman > 0.1)

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
mean(filt_en_cau_2_mets$spearman)
filt_en_cau_2_mets <- subset(filt_en_cau_2_mets, spearman > 0.1)

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
mean(filt_en_his_2_mets$spearman)
filt_en_his_2_mets <- subset(filt_en_his_2_mets, spearman > 0.1)

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
mean(filt_en_afa_2_mets$spearman)
filt_en_afa_2_mets <- subset(filt_en_afa_2_mets, spearman > 0.1)

#
#Make violin plot for all 3 population mesa to mets with EN without doing R > 0.1
en_afa_2mets <- data.frame(spearman=filt_en_afa_2_mets$spearman, 
                           prediction=rep("AFA", length(filt_en_afa_2_mets$spearman)))

en_his_2mets <- data.frame(spearman=filt_en_his_2_mets$spearman, 
                           prediction=rep("HIS", length(filt_en_his_2_mets$spearman)))

en_cau_2mets <- data.frame(spearman=filt_en_cau_2_mets$spearman, 
                           prediction=rep("CAU", length(filt_en_cau_2_mets$spearman)))

en_mesa_2_mets <- rbind(en_afa_2mets, en_his_2mets, en_cau_2mets)

library(ggplot2)
ggplot(en_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.2, color="black", lwd=1.2) +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.2), limits=c(-1.0, 1.0)) +
  theme(panel.grid.minor = element_blank()) + xlab("Elastic Net")

#
#Make violin plot for all 3 population mesa to mets with EN without doing R > 0.1
elnet_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
elnet_afa_2_mets$gene <- as.character(elnet_afa_2_mets$gene)
en_afa_2mets <- data.frame(spearman=elnet_afa_2_mets$spearman, 
                           prediction=rep("AFA", length(elnet_afa_2_mets$spearman)))

elnet_his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
elnet_his_2_mets$gene <- as.character(elnet_his_2_mets$gene)
en_his_2mets <- data.frame(spearman=elnet_his_2_mets$spearman, 
                           prediction=rep("HIS", length(elnet_his_2_mets$spearman)))

elnet_cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
elnet_cau_2_mets$gene <- as.character(elnet_cau_2_mets$gene)
en_cau_2mets <- data.frame(spearman=elnet_cau_2_mets$spearman, 
                           prediction=rep("CAU", length(elnet_cau_2_mets$spearman)))

en_mesa_2_mets <- rbind(en_afa_2mets, en_his_2mets, en_cau_2mets)
# Function to produce summary statistics (mean and +/- sd)

library(ggplot2)
ggplot(en_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.2, color="black", lwd=1.2) +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.2), limits=c(-1.0, 1.0)) +
  theme(panel.grid.minor = element_blank()) + xlab("Elastic Net")



#
##
###
##Make violin plot for all 3 population mesa to mets with ML without doing R > 0.1
library(dplyr)
#RF
#AFA
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
names(filt_rf_afa_2_mets) <- c("gene", "spearman")
mean(filt_rf_afa_2_mets$spearman)
filt_rf_afa_2_mets <- subset(filt_rf_afa_2_mets, spearman >= 0.1)


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
mean(filt_rf_his_2_mets$spearman)
filt_rf_his_2_mets <- subset(filt_rf_his_2_mets, spearman >= 0.1)

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
mean(filt_rf_cau_2_mets$spearman)
filt_rf_cau_2_mets <- subset(filt_rf_cau_2_mets, spearman >= 0.1)


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
mean(filt_rf_afhi_2_mets$spearman)
filt_rf_afhi_2_mets <- subset(filt_rf_afhi_2_mets, spearman >= 0.1)



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
mean(filt_rf_all_2_mets$spearman)
filt_rf_all_2_mets <- subset(filt_rf_all_2_mets, spearman >= 0.1)



#
#Make violin plot for all 3 population mesa to mets with RF without doing R > 0.1
rf_afa_2mets <- data.frame(spearman=filt_rf_afa_2_mets$spearman, 
                           prediction=rep("AFA", length(filt_rf_afa_2_mets$spearman)))

rf_his_2mets <- data.frame(spearman=filt_rf_his_2_mets$spearman, 
                           prediction=rep("HIS", length(filt_rf_his_2_mets$spearman)))

rf_cau_2mets <- data.frame(spearman=filt_rf_cau_2_mets$spearman, 
                           prediction=rep("CAU", length(filt_rf_cau_2_mets$spearman)))

rf_mesa_2_mets <- rbind(rf_afa_2mets, rf_his_2mets, rf_cau_2mets)

library(ggplot2)
ggplot(rf_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.2, color="black", lwd=1.2) +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.2), limits=c(-1.0, 1.0)) +
  theme(panel.grid.minor = element_blank()) + xlab("Random Forest")




#SVR
library(dplyr)
#AFA
svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_svr_all_chrom.txt", header = T)
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
mean(filt_svr_afa_2_mets$spearman)

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
mean(filt_svr_his_2_mets$spearman)


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
mean(filt_svr_cau_2_mets$spearman)


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
mean(filt_svr_afhi_2_mets$spearman)


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
mean(filt_svr_all_2_mets$spearman)

#
#Make violin plot for all 3 population mesa to mets with RF without doing R > 0.1
svr_afa_2mets <- data.frame(spearman=filt_svr_afa_2_mets$spearman, 
                           prediction=rep("AFA", length(filt_svr_afa_2_mets$spearman)))

svr_his_2mets <- data.frame(spearman=filt_svr_his_2_mets$spearman, 
                           prediction=rep("HIS", length(filt_svr_his_2_mets$spearman)))

svr_cau_2mets <- data.frame(spearman=filt_svr_cau_2_mets$spearman, 
                           prediction=rep("CAU", length(filt_svr_cau_2_mets$spearman)))

svr_mesa_2_mets <- rbind(svr_afa_2mets, svr_his_2mets, svr_cau_2mets)

library(ggplot2)
ggplot(svr_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.2, color="black", lwd=1.2) +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.2), limits=c(-1.0, 1.0)) +
  theme(panel.grid.minor = element_blank()) + xlab("SVR")



#KNN
library(dplyr)
#AFA
knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_knn_all_chrom.txt", header = T)
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
mean(filt_knn_afa_2_mets$spearman)

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
mean(filt_knn_his_2_mets$spearman)


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
mean(filt_knn_cau_2_mets$spearman)


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
mean(filt_knn_afhi_2_mets$spearman)


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
mean(filt_knn_all_2_mets$spearman)


#
#Make violin plot for all 3 population mesa to mets with KNN without doing R > 0.1
knn_afa_2mets <- data.frame(spearman=filt_knn_afa_2_mets$spearman, 
                           prediction=rep("AFA", length(filt_knn_afa_2_mets$spearman)))

knn_his_2mets <- data.frame(spearman=filt_knn_his_2_mets$spearman, 
                           prediction=rep("HIS", length(filt_knn_his_2_mets$spearman)))

knn_cau_2mets <- data.frame(spearman=filt_knn_cau_2_mets$spearman, 
                           prediction=rep("CAU", length(filt_knn_cau_2_mets$spearman)))

knn_mesa_2_mets <- rbind(knn_afa_2mets, knn_his_2mets, knn_cau_2mets)

library(ggplot2)
ggplot(knn_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.2, color="black", lwd=1.2) +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.2), limits=c(-1.0, 1.0)) +
  theme(panel.grid.minor = element_blank()) + xlab("KNN")



#############make violin plot of all data in table 2, and also do the model comparison for all the pops and AFA and AFHI
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

#width=1000, height=750

#compare ml and en across pops

all_en_rf <- inner_join(filt_en_all_2_mets, filt_rf_all_2_mets, by =c("gene_id"="gene"))
all_en_rf <- subset(all_en_rf, spearman.x > 0.1 & spearman.y > 0.1)
all_en_rf <- all_en_rf[,c(3,4)]
names(all_en_rf) <- c("x","y")

p1 <- ggplot(all_en_rf, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

all_en_svr <- inner_join(filt_en_all_2_mets, filt_svr_all_2_mets, by =c("gene_id"="gene"))
all_en_svr <- subset(all_en_svr, spearman.x > 0.1 & spearman.y > 0.1)
all_en_svr <- all_en_svr[,c(3,4)]
names(all_en_svr) <- c("x","y")

p2 <- ggplot(all_en_svr, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

all_en_knn <- inner_join(filt_en_all_2_mets, filt_knn_all_2_mets, by =c("gene_id"="gene"))
all_en_knn <- subset(all_en_knn, spearman.x > 0.1 & spearman.y > 0.1)
all_en_knn <- all_en_knn[,c(3,4)]
names(all_en_knn) <- c("x","y")

p3 <- ggplot(all_en_knn, aes(x=x, y=y)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("K Nearest Neighbor") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=820, height=600)

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
  scale_x_continuous(name = "Spearman Correlation") +
  geom_density(size=1.5)+ theme_classic(20) # +


#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=3200, height=2200

#head to head comparison of EN vs ML on performance in METS using AFA
library(dplyr)
#EN
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
filt_en_afa_2_mets <- filt_en_afa_2_mets[,c(1,3)]
names(filt_en_afa_2_mets) <- c("gene", "spearman")
filt_en_afa_2_mets <- subset(filt_en_afa_2_mets, spearman>0.1) #filtr by R >0.1
mean(filt_en_afa_2_mets$spearman)

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
names(filt_rf_afa_2_mets) <- c("gene", "spearman")
filt_rf_afa_2_mets <- subset(filt_rf_afa_2_mets, spearman>0.1) #filtr by R >0.1
mean(filt_rf_afa_2_mets$spearman)

#SVR
svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_svr_all_chrom.txt", header = T)
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
filt_svr_afa_2_mets <- subset(filt_svr_afa_2_mets, spearman>0.1) #filtr by R >0.1
mean(filt_svr_afa_2_mets$spearman)

#KNN
knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_knn_all_chrom.txt", header = T)
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
filt_knn_afa_2_mets <- subset(filt_knn_afa_2_mets, spearman>0.1) #filtr by R >0.1
mean(filt_knn_afa_2_mets$spearman)


library(dplyr)
en_rf <- inner_join(filt_en_afa_2_mets, filt_rf_afa_2_mets, by = c("gene"="gene"))
en_rf <- en_rf[,c(2,3)]
names(en_rf) <- c("en","rf")
cor.test(en_rf$en, en_rf$rf)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(filt_en_afa_2_mets, filt_svr_afa_2_mets, by = c("gene"="gene"))
en_svr <- en_svr[,c(2,3)]
names(en_svr) <- c("en","svr")
cor.test(en_svr$en, en_svr$svr)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(filt_en_afa_2_mets, filt_knn_afa_2_mets, by = c("gene"="gene"))
en_knn <- en_knn[,c(2,3)]
names(en_knn) <- c("en","knn")
cor.test(en_knn$en, en_knn$knn)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(0,1) + ylim(0,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)


#check and find genes in EN and RF AFA2METS models
en_rfboth <- inner_join(filt_en_afa_2_mets, filt_rf_afa_2_mets, by = c("gene"="gene"))
names(en_rfboth) <- c("gene", "en_spearman", "rf_spearman")
enplusrf_rf <- data.frame(spearman=en_rfboth$rf_spearman,prediction=rep("EN+RF RF Only",length(en_rfboth$rf_spearman)))
enplusrf_en <- data.frame(spearman=en_rfboth$en_spearman,prediction=rep("EN+RF EN Only",length(en_rfboth$en_spearman)))

rfonly <- anti_join(filt_rf_afa_2_mets, filt_en_afa_2_mets, by = c("gene" = "gene"))
rfonly <- anti_join(rfonly, en_rfboth, by = c("gene" = "gene")) #2nd check
rfonly_den <- data.frame(spearman=rfonly$spearman, prediction=rep("RF Only", length(rfonly$spearman)))

enonly <- anti_join(filt_en_afa_2_mets, filt_rf_afa_2_mets, by = c("gene" = "gene")) 
enonly <- anti_join(enonly, en_rfboth, by = c("gene" = "gene")) #2nd check
enonly_den <- data.frame(spearman=enonly$spearman, prediction=rep("EN Only", length(enonly$spearman)))

den_afa_2_mets <- rbind(rfonly_den, enonly_den, enplusrf_rf, enplusrf_en)

#find the mean of each ML group
library(plyr)
mu <- ddply(den_afa_2_mets, "prediction", summarise, grp.median=median(spearman))

# lwd = line thickness
ggplot(den_afa_2_mets, aes(x = spearman, color=prediction, lwd=1)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") +
  geom_density()+ theme_classic(20) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=prediction),linetype="longdash", lwd=1) +
  scale_color_manual(values = c("red","blue","orange","violet")) #+
#geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")