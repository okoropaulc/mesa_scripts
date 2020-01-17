#check the direction of effect for the hit genes in rf(TRPM4=ENSG00000130529.11) and en and svr(MCM3AP=ENSG00000160294.6)
#knn (TMEM50B=ENSG00000142188.12)
library(data.table)
#check gene name
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
# rf(ZFP36L1=ENSG00000185650.8, sox15=ENSG00000129194.3). and en (DESI1=ENSG00000100418.7)
#read in the prediction from both algs

pxcan <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/en_cau_rankplt5_predicted_expression.txt", header=T)
predrf <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_expr.txt",header=T)
predsvr <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_svr_pred_expr.txt",header=T)
predknn <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_knn_pred_expr.txt",header=T)


#read in the phenotype
thrombomodulin <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)
library(tidyverse)

#drop NA
#keep ID and rankplt5
thrombomodulin <- thrombomodulin[,c(1,5)]
thrombomodulin <- drop_na(thrombomodulin) #drop NA
thrombomodulin$sidno <- as.character(thrombomodulin$sidno)
colnames(thrombomodulin)[2] <- "phenotype" #rank_plt5

#remove one of the ID
pxcan$FID <- NULL

pxcan$IID <- as.character(pxcan$IID)
predrf$sampleidrf <-as.character(predrf$sampleidrf)
predsvr$sampleidsvr <- as.character(predsvr$sampleidsvr)
predknn$sampleidknn <- as.character(predknn$sampleidknn)

library(dplyr)

#merge pheno and pred expr on their common sample

pxcanmerged <- inner_join(thrombomodulin, pxcan, by = c("sidno" = "IID"))
rfmerged <- inner_join(thrombomodulin, predrf, by = c("sidno" = "sampleidrf"))
svrmerged <- inner_join(thrombomodulin, predsvr, by = c("sidno" = "sampleidsvr"))
knnmerged <- inner_join(thrombomodulin, predknn, by = c("sidno" = "sampleidknn"))

#pxcanmerged$sidno <- NULL

#extract the predictd expression for the genes
en_mcm3ap <- pxcanmerged[["ENSG00000160294.6"]]
en_desi1 <- pxcanmerged[["ENSG00000100418.7"]]
rf_trpm4 <-rfmerged[["ENSG00000130529.11"]]
rf_zfp36l1 <- rfmerged[["ENSG00000185650.8"]]
rf_sox15 <- rfmerged[["ENSG00000129194.3"]]
svr_mcm3ap <- svrmerged[["ENSG00000160294"]] #mcm3ap gene id does not have . in svr
knn_tmem50b <- knnmerged[["ENSG00000142188"]] #tmem50b gene id does not have . in knn
pheno <- rfmerged[["phenotype"]] #the phenotype is same for all merged df

#make a df and store pheno and genes of interest expression values. it will be easy to plot them that way
pheno_g <- data.frame(pheno=pheno,en_mcm3ap=en_mcm3ap,en_desi1=en_desi1,rf_trpm4=rf_trpm4,rf_zfp36l1=rf_zfp36l1,rf_sox15=rf_sox15,
                      svr_mcm3ap=svr_mcm3ap, knn_tmem50b=knn_tmem50b)

#plot them
#MCM3AP en
library(ggplot2)
ggplot(pheno_g, aes(x=en_mcm3ap, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") + xlim(-1.15,0.15) +
  theme_classic(20) + ggtitle("Elastic Net Gene MCM3AP") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
  #+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#MCM3AP svr
library(ggplot2)
ggplot(pheno_g, aes(x=svr_mcm3ap, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") + xlim(-0.85,0.5) +
  theme_classic(20) + ggtitle("SVR Gene MCM3AP") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#DESI1
library(ggplot2)
ggplot(pheno_g, aes(x=en_desi1, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") + xlim(-1.8,0.07) +
  theme_classic(20) + ggtitle("Elastic Net Gene DESI1") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#TRPM4
library(ggplot2)
ggplot(pheno_g, aes(x=rf_trpm4, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") + xlim(-2.0,2.2) +
  theme_classic(20) + ggtitle("Random Forest Gene TRPM4") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#ZFP36L1
library(ggplot2)
ggplot(pheno_g, aes(x=rf_zfp36l1, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") + xlim(-1.8,2.2) +
  theme_classic(20) + ggtitle("Random Forest Gene ZFP36L1") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#SOX15
library(ggplot2)
ggplot(pheno_g, aes(x=rf_sox15, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") + xlim(-2.0,5.0) +
  theme_classic(20) + ggtitle("Random Forest Gene SOX15") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#TMEM50B
library(ggplot2)
ggplot(pheno_g, aes(x=knn_tmem50b, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") + xlim(-3.1,4.1) +
  theme_classic(20) + ggtitle("KNN Gene TMEM50B") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")


# plot(rf_trpm4,pheno)
# plot(rf_zfp36l1, pheno)
# plot(rf_sox15, pheno)
# plot(en_mcm3ap, pheno)
# plot(en_desi1, pheno)

#compare the t statistic of rf and en on their overlap genes
# read in the filtered assoc files and use gene_name to overlap them
rfassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/rf_assoc_filtered_by_r2_0.01.txt", header=T)
rfassoc_cvr0.01$Gene_Name <- as.character(rfassoc_cvr0.01$Gene_Name)
enassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/en_assoc_filtered_by_r2_0.01.txt", header=T)
enassoc_cvr0.01$gene_name <- as.character(enassoc_cvr0.01$gene_name)

rfen <- inner_join(rfassoc_cvr0.01, enassoc_cvr0.01, by =c("Gene_Name"="gene_name"))
rfen <- rfen[,c(5,11)] #take out their t statistic
names(rfen) <- c("rf","en")
# #first remove decimals from the gene id
# for (i in 1:length(rf_afa$Gene_ID)){
#   rf_afa$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_afa$Gene_ID[i])
# } #just to remove the decimal places in the gene_id

library(ggplot2)
ggplot(rfen, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("Random Forest") + xlim(-6,6) + ylim(-6,6) +
  theme_classic(20) + ggtitle("t-statistic comparison") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)

#use ggpubr
library("ggpubr")
ggscatter(rfen, x = "en", y = "rf", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "Random Forest", 
          title = "t-statistic comparison",
          xlim = c(-6, 6), ylim = c(-6, 6)) + geom_abline(intercept = 0, slope = 1, color="blue")

