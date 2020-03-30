library(dplyr)
#AFA
en_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
en_afa$gene_id <- as.character(en_afa$gene_id)
en_afa <- en_afa[,c(1,10)]
names(en_afa) <- c("gene", "en")
en_afa <- subset(en_afa, en > -1)

#best_grid_rf_all_chrom.txt
rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_rf_all_chrom.txt", header = T)
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- rf_afa[,c(1,3)]
names(rf_afa) <- c("gene", "rf")
rf_afa <- subset(rf_afa, rf > -1)

svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header = T)
svr_afa$Gene_ID <- as.character(svr_afa$Gene_ID)
svr_afa <- svr_afa[,c(1,3)]
names(svr_afa) <- c("gene", "svr")
svr_afa <- subset(svr_afa, svr > -1)

knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header = T)
knn_afa$Gene_ID <- as.character(knn_afa$Gene_ID)
knn_afa <- knn_afa[,c(1,3)]
names(knn_afa) <- c("gene", "knn")
knn_afa <- subset(knn_afa, knn > -1)

# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_afa, rf_afa, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p1 <- ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_afa, svr_afa, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p2 <- ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_afa, knn_afa, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p3 <- ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(70)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)


####################
#plot violin plot of the cv r2 of en, rf, svr, knn
en <- data.frame(alg=rep("EN", nrow(en_afa)), r2=en_afa$en)
rf <- data.frame(alg=rep("RF", nrow(rf_afa)), r2=rf_afa$rf)
svr <- data.frame(alg=rep("SVR", nrow(svr_afa)), r2=svr_afa$svr)
knn <- data.frame(alg=rep("KNN", nrow(knn_afa)), r2=knn_afa$knn)

r2_comp <- rbind(en,rf,svr,knn)
# p4 <- ggplot(r2_comp, aes(x=alg, y=r2, color=alg, fill=alg)) + 
#   geom_violin(trim=F) + theme_classic(50) + xlab("") + theme(legend.position = "none") +
#   geom_boxplot(width=0.15,color="black") + ylab("CV Performance") + xlab("Model")

p4 <- ggplot(r2_comp, aes(x=alg, y=r2, fill=alg)) + geom_boxplot(lwd=0.5) + theme_light(70, base_line_size = 70/80) +
  xlab("Model") + ylab("CV Performance") + theme(legend.position = "none") + labs(title="D") +
  scale_y_continuous(breaks=seq(-1, 1.0, 0.25), limits=c(-1.0, 1.0)) #, fill=alg

#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=4000, height=3400




#HIS
en_his <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_HIS_model_summaries.txt", header = TRUE)
en_his$gene_id <- as.character(en_his$gene_id)
library(tidyverse)
en_his <- drop_na(en_his)
en_his <- subset(en_his, en_his$cv_R2_avg > -1)
en_his <- en_his[,c(1,10)]
names(en_his) <- c("gene","en")

rf_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_all_chrom.txt", header = T)
rf_his$Gene_ID <- as.character(rf_his$Gene_ID)
rf_his <- subset(rf_his, rf_his$CV_R2 > -1)
rf_his <- rf_his[,c(1,3)]
names(rf_his) <- c("gene", "rf")

svr_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_svr_all_chrom.txt", header = T)
svr_his$Gene_ID <- as.character(svr_his$Gene_ID)
svr_his <- subset(svr_his, svr_his$CV_R2 > -1)
svr_his <- svr_his[,c(1,3)]
names(svr_his) <- c("gene","svr")

knn_his <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_knn_all_chrom.txt", header = T)
knn_his$Gene_ID <- as.character(knn_his$Gene_ID)
knn_his <- subset(knn_his, knn_his$CV_R2 > -1)
knn_his <- knn_his[,c(1,3)]
names(knn_his) <- c("gene", "knn")

# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_his, rf_his, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p1 <- ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.2,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70) + labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_his, svr_his, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p2 <- ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.2,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70) + labs(title="B")# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_his, knn_his, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p3 <- ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.2,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(70) + labs(title="C")# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

####################
#plot violin plot of the cv r2 of en, rf, svr, knn
en <- data.frame(alg=rep("EN", nrow(en_his)), r2=en_his$en)
rf <- data.frame(alg=rep("RF", nrow(rf_his)), r2=rf_his$rf)
svr <- data.frame(alg=rep("SVR", nrow(svr_his)), r2=svr_his$svr)
knn <- data.frame(alg=rep("KNN", nrow(knn_his)), r2=knn_his$knn)

r2_comp <- rbind(en,rf,svr,knn)
# p4 <- ggplot(r2_comp, aes(x=alg, y=r2, color=alg, fill=alg)) + 
#   geom_violin(trim=F) + theme_classic(50) + xlab("") + theme(legend.position = "none") +
#   geom_boxplot(width=0.15,color="black") + ylab("CV Performance") + xlab("Model")

p4 <- ggplot(r2_comp, aes(x=alg, y=r2, fill=alg)) + geom_boxplot(lwd=0.5) + theme_light(70, base_line_size = 70/80) +
  xlab("Model") + ylab("CV Performance") + theme(legend.position = "none") + labs(title="D") +
  scale_y_continuous(breaks=seq(-1, 1.0, 0.25), limits=c(-0.75, 1.0)) #, fill=alg

#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=3500, height=2500


#CAU
en_cau <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_CAU_model_summaries.txt", header = TRUE)
en_cau$gene_id <- as.character(en_cau$gene_id)
library(tidyverse)
en_cau <- drop_na(en_cau)
#check for NA's in the avg_cv_R2 and drop it
en_cau <- subset(en_cau, en_cau$cv_R2_avg > -1)
en_cau <- en_cau[,c(1,10)]
names(en_cau) <- c("gene","en")

rf_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_rf_all_chrom.txt", header = T)
rf_cau$Gene_ID <- as.character(rf_cau$Gene_ID)
rf_cau <- subset(rf_cau, rf_cau$CV_R2 > -1)
rf_cau <- rf_cau[,c(1,3)]
names(rf_cau) <- c("gene", "rf")

svr_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_svr_all_chrom.txt", header = T)
svr_cau$Gene_ID <- as.character(svr_cau$Gene_ID)
svr_cau <- subset(svr_cau, svr_cau$CV_R2 > -1)
svr_cau <- svr_cau[,c(1,3)]
names(svr_cau) <- c("gene","svr")

knn_cau <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_knn_all_chrom.txt", header = T)
knn_cau$Gene_ID <- as.character(knn_cau$Gene_ID)
knn_cau <- subset(knn_cau, knn_cau$CV_R2 > -1)
knn_cau <- knn_cau[,c(1,3)]
names(knn_cau) <- c("gene", "knn")

# ML R2 comparisons against EN
library(dplyr)
en_rf <- inner_join(en_cau, rf_cau, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p1 <- ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(70) + labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_cau, svr_cau, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p2 <- ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(70) + labs(title="B")# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_cau, knn_cau, by = c("gene"="gene"))
library(ggplot2) #i used geom smooth to add the red regression line
p3 <- ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-0.1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(70) + labs(title="C")# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)


####################
#plot violin plot of the cv r2 of en, rf, svr, knn
en <- data.frame(alg=rep("EN", nrow(en_cau)), r2=en_cau$en)
rf <- data.frame(alg=rep("RF", nrow(rf_cau)), r2=rf_cau$rf)
svr <- data.frame(alg=rep("SVR", nrow(svr_cau)), r2=svr_cau$svr)
knn <- data.frame(alg=rep("KNN", nrow(knn_cau)), r2=knn_cau$knn)

r2_comp <- rbind(en,rf,svr,knn)
# p4 <- ggplot(r2_comp, aes(x=alg, y=r2, color=alg, fill=alg)) + 
#   geom_violin(trim=F) + theme_classic(50) + xlab("") + theme(legend.position = "none") +
#   geom_boxplot(width=0.15,color="black") + ylab("CV Performance") + xlab("Model")

p4 <- ggplot(r2_comp, aes(x=alg, y=r2, fill=alg)) + geom_boxplot(lwd=0.5) + theme_light(70, base_line_size = 70/80) +
  xlab("Model") + ylab("CV Performance") + theme(legend.position = "none") + labs(title="D") +
  scale_y_continuous(breaks=seq(-1, 1.0, 0.25), limits=c(-1.0, 1.0)) #, fill=alg


#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=3500, height=2500

#####AFHI
en_afhi <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFHI_model_summaries.txt", header=T)#it has no header
en_afhi <- en_afhi[,c(1,2,10)]
en_afhi$cv_R2_avg <- as.numeric(en_afhi$cv_R2_avg)
en_afhi$gene_id <- as.character(en_afhi$gene_id)


rf_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_rf_all_chrom.txt", header=T)
rf_afhi$CV_R2 <- as.numeric(rf_afhi$CV_R2)
rf_afhi$Gene_ID <- as.character(rf_afhi$Gene_ID)


svr_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_svr_all_chrom.txt", header=T)
svr_afhi$CV_R2 <- as.numeric(svr_afhi$CV_R2)
svr_afhi$Gene_ID <- as.character(svr_afhi$Gene_ID)


knn_afhi <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFHI_best_grid_knn_all_chrom.txt", header=T)
knn_afhi$CV_R2 <- as.numeric(knn_afhi$CV_R2)
knn_afhi$Gene_ID <- as.character(knn_afhi$Gene_ID)

# ML R2 comparisons against EN for AFHI
library(dplyr)
en_rf <- inner_join(en_afhi, rf_afhi, by = c("gene_id"="Gene_ID"))
en_rf <- en_rf[,c(3,5)]
names(en_rf) <- c("en", "rf")
library(ggplot2) #i used geom smooth to add the red regression line
p1 <- ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(50)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_svr <- inner_join(en_afhi, svr_afhi, by = c("gene_id"="Gene_ID"))
en_svr <- en_svr[,c(3,5)]
names(en_svr) <- c("en", "svr")
library(ggplot2) #i used geom smooth to add the red regression line
p2 <- ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("Support Vector") +
  theme_classic(50)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

en_knn <- inner_join(en_afhi, knn_afhi, by = c("gene_id"="Gene_ID"))
en_knn <- en_knn[,c(3,5)]
names(en_knn) <- c("en","knn")
library(ggplot2) #i used geom smooth to add the red regression line
p3 <- ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlim(-1,1) + ylim(-1,1) + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +
  theme_classic(50)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600)

#plot violin plot of the cv r2 of en, rf, svr, knn
en <- data.frame(alg=rep("EN", nrow(en_afhi)), r2=en_afhi$en)
rf <- data.frame(alg=rep("RF", nrow(rf_afhi)), r2=rf_afhi$rf)
svr <- data.frame(alg=rep("SVR", nrow(svr_afhi)), r2=svr_afhi$svr)
knn <- data.frame(alg=rep("KNN", nrow(knn_afhi)), r2=knn_afhi$knn)

r2_comp <- rbind(en,rf,svr,knn)
p4 <- ggplot(r2_comp, aes(x=alg, y=r2, color=alg, fill=alg)) + 
  geom_violin(trim=F) + theme_classic(50) + xlab("") + theme(legend.position = "none") +
  geom_boxplot(width=0.15,color="black") + ylab("CV Performance") + xlab("Model")

######################
###########               ALL

en_all <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", header=T)#it has no header
en_all <- en_all[,c(1,2,10)]
library(tidyverse)
en_all <- drop_na(en_all) #remove NA
en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)
en_all <- subset(en_all, cv_R2_avg > -0.5)
en_all$gene_id <- as.character(en_all$gene_id)


rf_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
rf_all$CV_R2 <- as.numeric(rf_all$CV_R2)
rf_all <- subset(rf_all, CV_R2 > -0.5)
rf_all$Gene_ID <- as.character(rf_all$Gene_ID)


svr_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
svr_all$CV_R2 <- as.numeric(svr_all$CV_R2)
svr_all <- subset(svr_all, CV_R2 > -0.5)
svr_all$Gene_ID <- as.character(svr_all$Gene_ID)


knn_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)
knn_all$CV_R2 <- as.numeric(knn_all$CV_R2)
knn_all <- subset(knn_all, CV_R2 > -0.5)
knn_all$Gene_ID <- as.character(knn_all$Gene_ID)



# ML R2 comparisons against EN for ALL
library(dplyr)
en_rf <- inner_join(en_all, rf_all, by = c("gene_id"="Gene_ID"))
en_rf <- en_rf[,c(3,5)]
names(en_rf) <- c("en", "rf")

library(ggplot2) #i used geom smooth to add the red regression line
p1 <- ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("Random Forest") + theme_classic(70) +
  xlim(-0.1,1) + ylim(-1,1) +
  labs(title="A")# + ggtitle("Cross Validation Perfromance") (save at width=900, height=600) + xlim(-1,1) + ylim(-1,1)

en_svr <- inner_join(en_all, svr_all, by = c("gene_id"="Gene_ID"))
en_svr <- en_svr[,c(3,5)]
names(en_svr) <- c("en", "svr")

library(ggplot2) #i used geom smooth to add the red regression line
p2 <- ggplot(en_svr, aes(x=en, y=svr)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("Support Vector") + theme_classic(70) + 
  xlim(-0.1,1) + ylim(-1,1) +
  labs(title="B")# before xlim(-0.1,1) + ylim(-1,1) + ggtitle("Cross Validation Perfromance") (save at width=900, height=600)+ xlim(-1,1) + ylim(-1,1)

en_knn <- inner_join(en_all, knn_all, by = c("gene_id"="Gene_ID"))
en_knn <- en_knn[,c(3,5)]
names(en_knn) <- c("en","knn")

library(ggplot2) #i used geom smooth to add the red regression line
p3 <- ggplot(en_knn, aes(x=en, y=knn)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("K-Nearest Neighbor") + theme_classic(70) +
  xlim(-0.1,1) + ylim(-1,1) +
  labs(title="C")# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600) + xlim(-1,1) + ylim(-1,1)


#####
#plot violin plot of the cv r2 of en, rf, svr, knn
en <- data.frame(alg=rep("EN", nrow(en_all)), r2=en_all$cv_R2_avg)
rf <- data.frame(alg=rep("RF", nrow(rf_all)), r2=rf_all$CV_R2)
svr <- data.frame(alg=rep("SVR", nrow(svr_all)), r2=svr_all$CV_R2)
knn <- data.frame(alg=rep("KNN", nrow(knn_all)), r2=knn_all$CV_R2)

r2_comp <- rbind(en,rf,svr,knn)
# p4 <- ggplot(r2_comp, aes(x=alg, y=r2, color=alg, fill=alg)) + 
#   geom_violin(trim=F) + theme_classic(70) + xlab("") + theme(legend.position = "none") +
#   geom_boxplot(width=0.15,color="black") + ylab("CV Performance") + xlab("Model") +
#   labs(title="D")


p4 <- ggplot(r2_comp, aes(x=alg, y=r2, fill=alg)) + geom_boxplot(lwd=0.5) + theme_light(70) +
  xlab("Model") + scale_y_continuous(breaks=seq(-0.5, 1.0, 0.25), limits=c(-0.5, 1.0)) +
  ylab("CV Performance") + theme(legend.position = "none") + labs(title="D")#, fill=alg

#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2) #save width=4000, height=3400
#grid.arrange(p1,p2,p3,nrow=1)
