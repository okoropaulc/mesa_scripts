"%&%" <- function(a,b) paste(a,b, sep = "")

for (chrom in c(7:14)){
  print(chrom)
  no <- as.character(chrom)
  knn_chr <- read.table("/home/paul/mesa_models/python_ml_models/results/knn_cv_chr" %&% no %&% ".txt", header = TRUE)
  png(filename = "/home/paul/mesa_models/python_ml_models/perf_plots/knn_chr" %&% no %&% ".png", width = 752, height = 519)
  hist(knn_chr$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "KNN Performance on MESA AFA chr" %&% no)
  dev.off()
  
  png(filename = "/home/paul/mesa_models/python_ml_models/perf_plots/rf_chr" %&% no %&% ".png", width = 752, height = 519)
  rf_chr <- read.table("/home/paul/mesa_models/python_ml_models/results/rf_cv_chr" %&% no %&% ".txt", header = TRUE)
  hist(rf_chr$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Random Forest Performance on MESA AFA chr" %&% no)
  dev.off()
  
  png(filename = "/home/paul/mesa_models/python_ml_models/perf_plots/svr_rbf_chr" %&% no %&% ".png", width = 752, height = 519)
  svr_rbf_chr <- read.table("/home/paul/mesa_models/python_ml_models/results/svr_rbf_cv_chr" %&% no %&% ".txt", header = TRUE)
  hist(svr_rbf_chr$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "SVR RBF Performance on MESA AFA chr" %&% no)
  dev.off()
  
  png(filename = "/home/paul/mesa_models/python_ml_models/perf_plots/svr_linear_chr" %&% no %&% ".png", width = 752, height = 519)
  svr_linear_chr <- read.table("/home/paul/mesa_models/python_ml_models/results/svr_linear_cv_chr" %&% no %&% ".txt", header = TRUE)
  hist(svr_linear_chr$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "SVR Linear Performance on MESA AFA chr" %&% no)
  dev.off()
  
  #Now check elnet chr22 R2 perf
  png(filename = "/home/paul/mesa_models/python_ml_models/perf_plots/elnet_chr" %&% no %&% ".png", width = 752, height = 519)
  elnet_chr <- read.table("/home/paul/mesa_models/split_mesa/results/full_AFA_chr" %&% no %&% "_model_summaries.txt", header = TRUE)
  hist(elnet_chr$cv_R2_avg, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Baseline(Elnet) Performance on MESA AFA chr" %&% no)
  dev.off()
}

#join all chromosomes to create one file for model results per ML algorithm
knn_model <- NULL
rf_model <- NULL
svr_rbf_model <- NULL
svr_linear_model <- NULL

for (chrom in 1:22) {
  no <- as.character(chrom)
  knn_model <- rbind(knn_model, read.table(file = "Z:/mesa_models/python_ml_models/results/knn_cv_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F))
  rf_model <- rbind(rf_model, read.table(file = "Z:/mesa_models/python_ml_models/results/rf_cv_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F))
  svr_rbf_model <- rbind(svr_rbf_model, read.table(file = "Z:/mesa_models/python_ml_models/results/svr_rbf_cv_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F))
  svr_linear_model <- rbind(svr_linear_model, read.table(file = "Z:/mesa_models/python_ml_models/results/svr_linear_cv_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F))
}

rf <- unique(rf_model)
knn <- unique(knn_model)

#write out the merged full chromosomes
write.table(knn, file="Z:/mesa_models/python_ml_models/results/knn_cv_full_AFA_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rf, file="Z:/mesa_models/python_ml_models/results/rf_cv_full_AFA_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr_rbf_model, file="Z:/mesa_models/python_ml_models/results/svr_rbf_cv_full_AFA_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr_linear_model, file="Z:/mesa_models/python_ml_models/results/svr_linear_cv_full_AFA_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#Compare python model
gencodev18 <- read.table(file = "/home/paul/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
#chromosome 22
knn_chr15 <- read.table("/home/paul/mesa_models/python_ml_models/results/knn_cv_chr15.txt", header = TRUE)
hist(knn_chr15$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "KNN Performance on MESA AFA chr15")

rf_chr15 <- read.table("/home/paul/mesa_models/python_ml_models/results/rf_cv_chr15.txt", header = TRUE)
hist(rf_chr15$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Random Forest Performance on MESA AFA chr15")

svr_rbf_chr15 <- read.table("/home/paul/mesa_models/python_ml_models/results/svr_rbf_cv_chr15.txt", header = TRUE)
hist(svr_rbf_chr15$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "SVR RBF Performance on MESA AFA chr15")

svr_linear_chr15 <- read.table("/home/paul/mesa_models/python_ml_models/results/svr_linear_cv_chr15.txt", header = TRUE)
hist(svr_linear_chr15$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "SVR Linear Performance on MESA AFA chr15")

#Now check elnet chr22 R2 perf
elnet_chr15 <- read.table("/home/paul/mesa_models/split_mesa/results/full_AFA_chr15_model_summaries.txt", header = TRUE)
hist(elnet_chr15$cv_R2_avg, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Baseline(Elnet) Performance on MESA AFA chr15")

#full elnet AFA chrom
elnet_all <- read.table(file = "Z:/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)

#plot the distribution of the ML models CV R2
hist(elnet_all$cv_R2_avg, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Baseline (Elnet) Performance on MESA AFA all chr")
hist(svr_linear_model$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "SVR Linear Performance on MESA AFA all chr")
hist(svr_rbf_model$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "SVR RBF Performance on MESA AFA all chr")
hist(rf_model$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Random Forest Performance on MESA AFA all chr")
hist(knn_model$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "KNN Performance on MESA AFA all chr")

#read in saved files
rf <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/rf_cv_full_AFA_chr.txt", header = TRUE)
elnet <- read.table(file = "/home/paul/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)

#sort by gene_id
rf$Gene_ID <- as.character(rf$Gene_ID)
elnet$gene_id <- as.character(elnet$gene_id)

#sort the dataframe by gene id
rf <- rf[order(rf$Gene_ID),]
elnet <- elnet[order(elnet$gene_id),]

el_rf <- cbind(elnet[,c(1,10)], rf[,1:2])
library(tidyverse)
ggplot(el_rf, aes(x=cv_R2_avg,y=CV_R2)) + ggtitle("Elastic Net vs Random Forest") + ylab("RF R2") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue')


#rf_model$Gene_ID <- as.character(rf_model$Gene_ID)
#el_rf <- inner_join(elnet_all, rf_model, by= c("gene_id" = "Gene_ID"))

#before filtering, they both have 9623 genes

#filter the R2 to be within -1 and 1
rf_filt <- subset(rf, CV_R2 >= -1)
elnet_filt <- subset(elnet, cv_R2_avg >= -1)
# after filtering, random forest has 9622, while elnet has 9609 genes left

rf_filt <- rf_filt[,1:2]
elnet_filt <- elnet_filt[,c(1,10)]

library(dplyr)

el_rf_filt <- inner_join(elnet_filt, rf_filt, by = c("gene_id" = "Gene_ID"))
names(el_rf_filt) <- c("gene_id", "elnet_cv_R2", "rf_cv_R2")

ggplot(el_rf_filt, aes(x=elnet_cv_R2, y=rf_cv_R2)) + ggtitle("Elastic Net vs Random Forest (outliers (< -1) filtered off)") + 
  ylab("RF R2") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

#check the correlation of the two: elnet_filt and rf_filt

library("ggpubr")
ggscatter(el_rf_filt, x = "elnet_cv_R2", y = "rf_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net R2", ylab = "RF R2", title = "Elastic Net vs Random Forest (outliers (< -1) filtered off)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#KNN
knn <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/knn_cv_full_AFA_chr.txt", header = TRUE)
knn$Gene_ID <- as.character(knn$Gene_ID)
knn$time.s. <- NULL

el_knn_filt <- inner_join(elnet_filt, knn, by = c("gene_id" = "Gene_ID"))
names(el_knn_filt) <- c("gene_id", "elnet_cv_R2", "knn_cv_R2")

ggplot(el_knn_filt, aes(x=elnet_cv_R2, y=knn_cv_R2)) + ggtitle("Elastic Net vs KNN (outliers (< -1) filtered off)") + 
  ylab("KNN R2") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

ggscatter(el_knn_filt, x = "elnet_cv_R2", y = "knn_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net R2", ylab = "KNN R2", title = "Elastic Net vs KNN (outliers (< -1) filtered off)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR RBF
svr_rbf <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/svr_rbf_cv_full_AFA_chr.txt", header = TRUE)
svr_rbf$Gene_ID <- as.character(svr_rbf$Gene_ID)
svr_rbf$time.s. <- NULL

el_svr_rbf_filt <- inner_join(elnet_filt, svr_rbf, by = c("gene_id" = "Gene_ID"))
names(el_svr_rbf_filt) <- c("gene_id", "elnet_cv_R2", "svr_rbf_cv_R2")

ggplot(el_svr_rbf_filt, aes(x=elnet_cv_R2, y=svr_rbf_cv_R2)) + ggtitle("Elastic Net vs SVR RBF (outliers (< -1) filtered off)") + 
  ylab("SVR RBF R2") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

ggscatter(el_svr_rbf_filt, x = "elnet_cv_R2", y = "svr_rbf_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net R2", ylab = "SVR RBF R2", title = "Elastic Net vs SVR RBF (outliers (< -1) filtered off)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR RBF
svr_linear <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/svr_linear_cv_full_AFA_chr.txt", header = TRUE)
svr_linear$Gene_ID <- as.character(svr_linear$Gene_ID)
svr_linear$time.s. <- NULL

el_svr_linear_filt <- inner_join(elnet_filt, svr_linear, by = c("gene_id" = "Gene_ID"))
names(el_svr_linear_filt) <- c("gene_id", "elnet_cv_R2", "svr_linear_cv_R2")

ggplot(el_svr_linear_filt, aes(x=elnet_cv_R2, y=svr_linear_cv_R2)) + ggtitle("Elastic Net vs SVR Linear (outliers (< -1) filtered off)") + 
  ylab("SVR Linear R2") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

ggscatter(el_svr_linear_filt, x = "elnet_cv_R2", y = "svr_linear_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net R2", ylab = "SVR Linear R2", title = "Elastic Net vs SVR Linear (outliers (< -1) filtered off)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


# Model performance on test data set

"%&%" <- function(a,b) paste(a,b, sep = "")

knn <- NULL
rf <- NULL
svr_rbf <- NULL
svr_linear <- NULL
#el <- NULL
pop <- "HIS"
i <- "METS"

for (chrom in 1:22) {
  no <- as.character(chrom)
  knn <- rbind(knn, read.table(file = "/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_knn_cor_test_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  rf <- rbind(rf, read.table(file = "/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_rf_cor_test_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  svr_rbf <- rbind(svr_rbf, read.table(file = "/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_svr_rbf_cor_test_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  svr_linear <- rbind(svr_linear, read.table(file = "/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_svr_linear_cor_test_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  
  #el <- rbind(el, read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_elnet_cor_test_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(knn, file="/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_knn_cor_test_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rf, file="/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_rf_cor_test_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr_rbf, file="/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_svr_rbf_cor_test_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr_linear, file="/home/paul/mesa_models/python_ml_models/results/" %&% pop %&% "_2_" %&% i %&% "_svr_linear_cor_test_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(el, file="/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_elnet_cor_test_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#checking pyelnet and rf
el$gene_id <- as.character(el$gene_id)
el <- el[,c(1,9)]
rfafa2his <- rf_afa2his[,c(1,8)]
elrf <- inner_join(el, rfafa2his, by = c("gene_id" = "gene_id"))
elrf <- elrf[,c(2,3)]
names(elrf) <- c("py_elnet", "rf")
ggplot(elrf, aes(x=py_elnet, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)") + 
  ylab("RF") + xlab("Python Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#check py rf and svr_linear
svrl <- svr_linear_afa2his[,c(1,8)]
rf_svrl <- inner_join(rfafa2his, svrl, by = c("gene_id" = "gene_id"))
rf_svrl <- rf_svrl[,c(2,3)]
names(rf_svrl) <- c("rf", "svrl")
ggplot(rf_svrl, aes(x=svrl, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)") + 
  ylab("RF") + xlab("SVR Linear") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#ML performances on test dataset
#HIS
#RF
#elnet_afa2his <- read.table(file = "/home/ryan/MESA_files/prediction_analysis/Spearman.MESA.AFA.on.MESAHIS.txt", header = T)
elnet_afa2his <- fread(file = "/home/ryan/MESA_files/prediction_analysis/redo/spearman.MESA_AFA_on_MESA_HIS.txt", header = T)
elnet_afa2his <- subset(elnet_afa2his, spearman != 0 | V3 != 0)
elnet_afa2his$V3 <- NULL
hist(elnet_afa2his$spearman, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) R Elastic Net")
elnet_afa2his$gene <- as.character(elnet_afa2his$gene)
rf_afa2his <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_HIS_rf_cor_test_all_chr.txt", header = T)
hist(rf_afa2his$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) RF")
rf_afa2his$gene_id <- as.character(rf_afa2his$gene_id)
for (i in 1:length(rf_afa2his$gene_id)){
  rf_afa2his$gene_id[i] <- gsub('\\.[0-9]+','',rf_afa2his$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_rf_filt <- inner_join(elnet_afa2his, rf_afa2his, by = c("gene" = "gene_id"))
el_rf_filt <- el_rf_filt[,c(2,9)]
names(el_rf_filt) <- c("elnet", "rf")
hist(el_rf_filt$rf, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) RF")

ggplot(el_rf_filt, aes(x=elnet, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)") + 
  ylab("RF") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_rf_filt, x = "elnet", y = "rf", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "RF", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#KNN
#elnet_afa2his <- read.table(file = "/home/ryan/MESA_files/prediction_analysis/Spearman.MESA.AFA.on.MESAHIS.txt", header = T)
elnet_afa2his <- fread(file = "/home/ryan/MESA_files/prediction_analysis/redo/spearman.MESA_AFA_on_MESA_HIS.txt", header = T)
elnet_afa2his <- subset(elnet_afa2his, spearman != 0 | V3 != 0)
elnet_afa2his$V3 <- NULL
elnet_afa2his$gene <- as.character(elnet_afa2his$gene)
knn_afa2his <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_HIS_knn_cor_test_all_chr.txt", header = T)
hist(knn_afa2his$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) KNN")
knn_afa2his$gene_id <- as.character(knn_afa2his$gene_id)
for (i in 1:length(knn_afa2his$gene_id)){
  knn_afa2his$gene_id[i] <- gsub('\\.[0-9]+','',knn_afa2his$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_knn_filt <- inner_join(elnet_afa2his, knn_afa2his, by = c("gene" = "gene_id"))
el_knn_filt <- el_knn_filt[,c(2,9)]
names(el_knn_filt) <- c("elnet", "knn")
hist(el_knn_filt$knn, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) KNN")


ggplot(el_knn_filt, aes(x=elnet, y=knn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_knn_filt, x = "elnet", y = "knn", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR_RBF
svr_rbf_afa2his <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_HIS_svr_rbf_cor_test_all_chr.txt", header = T)
hist(svr_rbf_afa2his$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) SVR RBF")
svr_rbf_afa2his$gene_id <- as.character(svr_rbf_afa2his$gene_id)
for (i in 1:length(svr_rbf_afa2his$gene_id)){
  svr_rbf_afa2his$gene_id[i] <- gsub('\\.[0-9]+','',svr_rbf_afa2his$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_svr_rbf_filt <- inner_join(elnet_afa2his, svr_rbf_afa2his, by = c("gene" = "gene_id"))
el_svr_rbf_filt <- el_svr_rbf_filt[,c(2,9)]
names(el_svr_rbf_filt) <- c("elnet", "svr_rbf")

hist(el_svr_rbf_filt$svr_rbf, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) SVR RBF")

ggplot(el_svr_rbf_filt, aes(x=elnet, y=svr_rbf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)") + 
  ylab("SVR RBF") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_svr_rbf_filt, x = "elnet", y = "svr_rbf", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "svr_rbf", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR_Linear
svr_linear_afa2his <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_HIS_svr_linear_cor_test_all_chr.txt", header = T)
hist(svr_linear_afa2his$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) SVR Linear")
svr_linear_afa2his$gene_id <- as.character(svr_linear_afa2his$gene_id)
for (i in 1:length(svr_linear_afa2his$gene_id)){
  svr_linear_afa2his$gene_id[i] <- gsub('\\.[0-9]+','',svr_linear_afa2his$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_svr_linear_filt <- inner_join(elnet_afa2his, svr_linear_afa2his, by = c("gene" = "gene_id"))
el_svr_linear_filt <- el_svr_linear_filt[,c(2,9)]
names(el_svr_linear_filt) <- c("elnet", "svr_linear")

hist(el_svr_linear_filt$svr_linear, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to HIS) SVR Linear")


ggplot(el_svr_linear_filt, aes(x=elnet, y=svr_linear)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)") + 
  ylab("SVR Linear") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_svr_linear_filt, x = "elnet", y = "svr_linear", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "svr_linear", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to HIS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")





#CAU
#RF
#elnet_afa2cau <- read.table(file = "/home/ryan/MESA_files/prediction_analysis/Spearman.MESA.AFA.on.MESACAU.txt", header = T)
elnet_afa2cau <- fread(file = "/home/ryan/MESA_files/prediction_analysis/redo/spearman.MESA_AFA_on_MESA_CAU.txt", header = T)
elnet_afa2cau <- subset(elnet_afa2cau, spearman != 0 | V3 != 0)
elnet_afa2cau$V3 <- NULL
hist(elnet_afa2cau$spearman, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) Elastic Net")
elnet_afa2cau$gene <- as.character(elnet_afa2cau$gene)

rf_afa2cau <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_CAU_rf_cor_test_all_chr.txt", header = T)
hist(rf_afa2cau$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) RF")
rf_afa2cau$gene_id <- as.character(rf_afa2cau$gene_id)
for (i in 1:length(rf_afa2cau$gene_id)){
  rf_afa2cau$gene_id[i] <- gsub('\\.[0-9]+','',rf_afa2cau$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_rf_filt <- inner_join(elnet_afa2cau, rf_afa2cau, by = c("gene" = "gene_id"))
el_rf_filt <- el_rf_filt[,c(2,9)]
names(el_rf_filt) <- c("elnet", "rf")
hist(el_rf_filt$rf, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) RF")


ggplot(el_rf_filt, aes(x=elnet, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)") + 
  ylab("RF") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_rf_filt, x = "elnet", y = "rf", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "RF", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#KNN
knn_afa2cau <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_CAU_knn_cor_test_all_chr.txt", header = T)
hist(knn_afa2cau$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) KNN")

knn_afa2cau$gene_id <- as.character(knn_afa2cau$gene_id)
for (i in 1:length(knn_afa2cau$gene_id)){
  knn_afa2cau$gene_id[i] <- gsub('\\.[0-9]+','',knn_afa2cau$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_knn_filt <- inner_join(elnet_afa2cau, knn_afa2cau, by = c("gene" = "gene_id"))
el_knn_filt <- el_knn_filt[,c(2,9)]
names(el_knn_filt) <- c("elnet", "knn")
hist(el_knn_filt$knn, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) KNN")


ggplot(el_knn_filt, aes(x=elnet, y=knn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_knn_filt, x = "elnet", y = "knn", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#SVR_RBF
svr_rbf_afa2cau <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_CAU_svr_rbf_cor_test_all_chr.txt", header = T)
hist(svr_rbf_afa2cau$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) SVR RBF")
svr_rbf_afa2cau$gene_id <- as.character(svr_rbf_afa2cau$gene_id)
for (i in 1:length(svr_rbf_afa2cau$gene_id)){
  svr_rbf_afa2cau$gene_id[i] <- gsub('\\.[0-9]+','',svr_rbf_afa2cau$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_svr_rbf_filt <- inner_join(elnet_afa2cau, svr_rbf_afa2cau, by = c("gene" = "gene_id"))
el_svr_rbf_filt <- el_svr_rbf_filt[,c(2,9)]
names(el_svr_rbf_filt) <- c("elnet", "svr_rbf")
hist(el_svr_rbf_filt$svr_rbf, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) SVR RBF")


ggplot(el_svr_rbf_filt, aes(x=elnet, y=svr_rbf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)") + 
  ylab("SVR RBF") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_svr_rbf_filt, x = "elnet", y = "svr_rbf", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "svr_rbf", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR_Linear
svr_linear_afa2cau <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_CAU_svr_linear_cor_test_all_chr.txt", header = T)
hist(svr_linear_afa2cau$spearman_yobs_vs_ypred..d., breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) SVR Linear")
svr_linear_afa2cau$gene_id <- as.character(svr_linear_afa2cau$gene_id)
for (i in 1:length(svr_linear_afa2cau$gene_id)){
  svr_linear_afa2cau$gene_id[i] <- gsub('\\.[0-9]+','',svr_linear_afa2cau$gene_id[i])
} #just to remove the decimal places in the gene_id

library(dplyr)
library(ggplot2)

el_svr_linear_filt <- inner_join(elnet_afa2cau, svr_linear_afa2cau, by = c("gene" = "gene_id"))
el_svr_linear_filt <- el_svr_linear_filt[,c(2,9)]
names(el_svr_linear_filt) <- c("elnet", "svr_linear")
hist(el_svr_linear_filt$svr_linear, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to CAU) SVR Linear")


ggplot(el_svr_linear_filt, aes(x=elnet, y=svr_linear)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)") + 
  ylab("SVR Linear") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

library("ggpubr")
ggscatter(el_svr_linear_filt, x = "elnet", y = "svr_linear", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "svr_linear", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to CAU)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")




#METS
library(data.table)
library(dplyr)
library(ggplot2)
#RF
el_afa2mets<-fread("/home/paul/spearman_AFA_2_METS.txt", header = T, sep = '\t') %>% select(gene,spearman)
hist(el_afa2mets$spearman, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (CAU to METS) R Elastic Net")
el_afa2mets$gene <- as.character(el_afa2mets$gene)

rf <- read.table(file="/home/paul/mesa_models/python_ml_models/results/HIS_2_METS_rf_cor_test_all_chr.txt", header = T, sep = "\t")

rf$gene_id <- as.character(rf$gene_id)
rf <- rf[,c(1,9)]
elrf <- inner_join(el_afa2mets, rf, by = c("gene" = "gene_id"))
names(elrf) <- c("gene","elnet","rf")
hist(elrf$elnet, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (HIS to METS) R Elastic Net")

ggplot(elrf, aes(x=elnet, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("RF") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

hist(elrf$rf, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (HIS to METS) RF")

library("ggpubr")
ggscatter(elrf, x = "elnet", y = "rf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "RF", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#KNN
el_afa2mets<-fread("/home/paul/spearman_AFA_2_METS.txt", header = T, sep = '\t') %>% select(gene,spearman)
el_afa2mets$gene <- as.character(el_afa2mets$gene)

knn <- read.table(file="/home/paul/mesa_models/python_ml_models/results/HIS_2_METS_knn_cor_test_all_chr.txt", header = T, sep = "\t")
knn$gene_id <- as.character(knn$gene_id)
knn <- knn[,c(1,9)]
elknn <- inner_join(el_afa2mets, knn, by = c("gene" = "gene_id"))
names(elknn) <- c("gene","elnet","knn")
ggplot(elknn, aes(x=elnet, y=knn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

hist(elknn$knn, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (HIS to METS) KNN")

ggscatter(elknn, x = "elnet", y = "knn", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR Linear
el_afa2mets<-fread("/home/paul/spearman_AFA.txt", header = T, sep = '\t') %>% select(gene,spearman)
el_afa2mets$gene <- as.character(el_afa2mets$gene)
svr_linear <- read.table(file="/home/paul/mesa_models/python_ml_models/results/HIS_2_METS_svr_linear_cor_test_all_chr.txt", header = T, sep = "\t")
svr_linear$gene_id <- as.character(svr_linear$gene_id)
svr_linear <- svr_linear[,c(1,9)]
elsvrl <- inner_join(el_afa2mets, svr_linear, by = c("gene" = "gene_id"))
names(elsvrl) <- c("gene","elnet","svrl")
ggplot(elsvrl, aes(x=elnet, y=svrl)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR Linear") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

hist(elsvrl$svrl, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (HIS to METS) SVR Linear")

ggscatter(elsvrl, x = "elnet", y = "svrl", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR Linear", title = "Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR RBF
el_afa2mets<-fread("/home/paul/spearman_AFA.txt", header = T, sep = '\t') %>% select(gene,spearman)
el_afa2mets$gene <- as.character(el_afa2mets$gene)

svr_rbf <- read.table(file="/home/paul/mesa_models/python_ml_models/results/HIS_2_METS_svr_rbf_cor_test_all_chr.txt", header = T, sep = "\t")
svr_rbf$gene_id <- as.character(svr_rbf$gene_id)
svr_rbf <- svr_rbf[,c(1,9)]
elsvr <- inner_join(el_afa2mets, svr_rbf, by = c("gene" = "gene_id"))
names(elsvr) <- c("gene","elnet","svr")
ggplot(elsvr, aes(x=elnet, y=svr)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR RBF") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

hist(elsvr$svr, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (HIS to METS) SVR RBF")
#hist(elsvr$elnet, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "Spearman", main = "Observed vs Predicted (AFA to METS) Elastic Net")

ggscatter(elsvr, x = "elnet", y = "svr", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR RBF", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS) to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#Merge all the pc adjusted gene expression for all chromosomes
pc <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/CAU_pc_adjusted_gene_expr_chr10.txt", header = T, sep = "\t")
pc2 <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/CAU_pc_adjusted_gene_expr_chr11.txt", header = T, sep = "\t")

#pcbind1 <- cbind(pc, pc2[2:length(pc2)])
pcbind <- cbind(pc, pc2[2:length(pc2)])
rn <- pc$X
pcbind$X <- NULL
colum <- names(pcbind)
pcbind <- as.data.frame(transpose(pcbind))
colnames(pcbind) <- rn
df <- as.data.frame(colum)
pcbind <- cbind(df, pcbind)

"%&%" <- function(a,b) paste(a,b, sep = "")
i <- "AFA"
pcadj <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/" %&% i %&% "_pc_adjusted_gene_expr_chr1.txt", header = T, stringsAsFactors = F, sep = "\t")

for (chrom in 2:22) {
  no <- as.character(chrom)
  pcx <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/" %&% i %&% "_pc_adjusted_gene_expr_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t")
  pcadj <- cbind(pcadj, pcx[2:length(pcx)])
}
id <- pcadj$X
pcadj$X <- NULL
gene_id <- names(pcadj)
pcadj <- as.data.frame(transpose(pcadj))
colnames(pcadj) <- id
df_col <- as.data.frame(gene_id)
pcadj <- cbind(df_col, pcadj)

#write out the merged full chromosomes
write.table(pcadj, file="/home/paul/mesa_models/python_ml_models/results/" %&% i %&% "_pc_adjusted_gene_expr_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)

pcs <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/CAU_pc_adjusted_gene_expr_all_chr.txt", header = T, sep = "\t")
pc_h <- fread(file = "/home/paul/mesa_models/python_ml_models/results/HIS_pc_adjusted_gene_expr_all_chr.txt", header = T, sep = "\t")
pc_a <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_pc_adjusted_gene_expr_all_chr.txt", header = T, sep = "\t",stringsAsFactors = F)

cau_gex <- read.table(file = "/home/paul/mesa_models/meqtl_sorted_CAU_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
cadj <- subset(pcs, pcs$gene_id == "ENSG00000116729.9")
cobs <- subset(cau_gex, cau_gex$PROBE_ID == "ENSG00000116729.9")
x <- as.numeric(cadj[1,])
y <- as.numeric(cobs[1,])
cor.test(x,y)
plot(x,y)


#Compare models performance on two pops eg.. elnet perf on afa2mets vs cau2mets
library(dplyr)
library(ggplot2)
library("ggpubr")

el_afa2mets<-fread("/home/paul/spearman_HIS_2_METS.txt", header = T, sep = '\t') %>% select(gene,spearman)
el_afa2mets$gene <- as.character(el_afa2mets$gene)
el_cau2mets<-fread("/home/paul/spearman_CAU_2_METS.txt", header = T, sep = '\t') %>% select(gene,spearman)
el_cau2mets$gene <- as.character(el_cau2mets$gene)
elel <- inner_join(el_afa2mets, el_cau2mets, by = c("gene" = "gene"))
names(elel) <- c("gene","afa2mets","cau2mets")


rf_afa2mets <- read.table(file="/home/paul/mesa_models/python_ml_models/results/AFA_2_METS_svr_linear_cor_test_all_chr.txt", header = T, sep = "\t")
rf_afa2mets <- rf_afa2mets[,c(1,9)]
rf_afa2mets$gene_id <- as.character(rf_afa2mets$gene_id)
rf_cau2mets <- read.table(file="/home/paul/mesa_models/python_ml_models/results/CAU_2_METS_svr_linear_cor_test_all_chr.txt", header = T, sep = "\t")
rf_cau2mets <- rf_cau2mets[,c(1,9)]
rf_cau2mets$gene_id <- as.character(rf_cau2mets$gene_id)
rfrf <- inner_join(rf_afa2mets, rf_cau2mets, by = c("gene_id" = "gene_id"))
names(rfrf) <- c("gene","afa2mets","cau2mets")

ggplot(rfrf, aes(x=afa2mets, y=cau2mets)) + ggtitle("Spearman Corr of AFA2METS vs CAU2METS (SVR Linear)") + 
  ylab("CAU2METS") + xlab("AFA2METS") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(rfrf, x = "afa2mets", y = "cau2mets", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "AFA2METS", ylab = "CAU2METS", title = "Spearman Corr of AFA2METS vs CAU2METS (SVR Linear)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")
