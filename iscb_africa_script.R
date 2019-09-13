#the svr used here is svr rbf
library(dplyr)

gencode <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode <- subset(gencode, gencode$gene_type == "protein_coding")
gencode <- gencode[,c(2,3)]
gencode$gene_name <- as.character(gencode$gene_name)
gencode2 <- gencode

#AFA2METS
afa_svr <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_svr_rbf_cor_test_all_chr.txt", header=T,sep="\t" )
afa_svr <- afa_svr[,c(1,9)]
afa_svr$gene_id <- as.character(afa_svr$gene_id)

afa_rf <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_rf_cor_test_all_chr.txt", header=T,sep="\t" )
afa_rf <- afa_rf[,c(1,9)]
afa_rf$gene_id <- as.character(afa_rf$gene_id)

afa_knn <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_knn_cor_test_all_chr.txt", header=T,sep="\t" )
afa_knn <- afa_knn[,c(1,9)]
afa_knn$gene_id <- as.character(afa_knn$gene_id)

ach_ar <- inner_join(ach, afa_rf, by = c('gene' = 'gene_id'))
ach_ars <- inner_join(ach_ar, afa_svr, by = c('gene' = 'gene_id'))
ach_arsk <- inner_join(ach_ars, afa_knn, by = c('gene' = 'gene_id'))

#AFA CV
afa_en_cv <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = T)
afa_en_cv <- afa_en_cv[,c(1,2,10)]
afa_en_cv$gene_name <- as.character(afa_en_cv$gene_name)
afa_en_cv <- subset(afa_en_cv, afa_en_cv$cv_R2_avg >= 0.7)

afa_svr_cv <- read.table(file="Z:/data/mesa_models/python_ml_models/results/svr_rbf_cv_full_AFA_chr.txt", header=T, sep="\t")
afa_svr_cv$Gene_ID <- as.character(afa_svr_cv$Gene_ID)
afa_svr_cv <- inner_join(gencode, afa_svr_cv, by = c("gene_id" = "Gene_ID"))
afa_svr_cv <-subset(afa_svr_cv, afa_svr_cv$CV_R2 >= 0.1)

afa_rf_cv <- read.table(file="Z:/data/mesa_models/python_ml_models/results/rf_cv_full_AFA_chr.txt", header=T, sep="\t")
afa_rf_cv$Gene_ID <- as.character(afa_rf_cv$Gene_ID)
afa_rf_cv <- inner_join(gencode, afa_rf_cv, by = c("gene_id" = "Gene_ID"))
afa_rf_cv <-subset(afa_rf_cv, afa_rf_cv$CV_R2 >= 0.7)
afa_rf_cv$gene_name <- as.character(afa_rf_cv$gene_name)

afa_knn_cv <- read.table(file="Z:/data/mesa_models/python_ml_models/results/knn_cv_full_AFA_chr.txt", header=T, sep="\t")
afa_knn_cv$Gene_ID <- as.character(afa_knn_cv$Gene_ID)
afa_knn_cv <- inner_join(gencode, afa_knn_cv, by = c("gene_id" = "Gene_ID"))
afa_knn_cv <-subset(afa_knn_cv, afa_knn_cv$CV_R2 >= 0.7)

#HIS2METS
his_svr <- read.table(file="Z:/data/mesa_models/python_ml_models/results/HIS_2_METS_svr_rbf_cor_test_all_chr.txt", header=T,sep="\t" )
his_svr <- his_svr[,c(1,9)]
his_svr$gene_id <- as.character(his_svr$gene_id)

his_rf <- read.table(file="Z:/data/mesa_models/python_ml_models/results/HIS_2_METS_rf_cor_test_all_chr.txt", header=T,sep="\t" )
his_rf <- his_rf[,c(1,9)]
his_rf$gene_id <- as.character(his_rf$gene_id)

his_knn <- read.table(file="Z:/data/mesa_models/python_ml_models/results/HIS_2_METS_knn_cor_test_all_chr.txt", header=T,sep="\t" )
his_knn <- his_knn[,c(1,9)]
his_knn$gene_id <- as.character(his_knn$gene_id)

ach_arsk_hr <- inner_join(ach_arsk, his_rf, by = c('gene' = 'gene_id'))
ach_arsk_hrs <- inner_join(ach_arsk_hr, his_svr, by = c('gene' = 'gene_id'))
ach_arsk_hrsk <- inner_join(ach_arsk_hrs, his_knn, by = c('gene' = 'gene_id'))

#CAU2METS
cau_svr <- read.table(file="Z:/data/mesa_models/python_ml_models/results/CAU_2_METS_svr_rbf_cor_test_all_chr.txt", header=T,sep="\t" )
cau_svr <- cau_svr[,c(1,9)]
cau_svr$gene_id <- as.character(cau_svr$gene_id)

cau_rf <- read.table(file="Z:/data/mesa_models/python_ml_models/results/CAU_2_METS_rf_cor_test_all_chr.txt", header=T,sep="\t" )
cau_rf <- cau_rf[,c(1,9)]
cau_rf$gene_id <- as.character(cau_rf$gene_id)

cau_knn <- read.table(file="Z:/data/mesa_models/python_ml_models/results/CAU_2_METS_knn_cor_test_all_chr.txt", header=T,sep="\t" )
cau_knn <- cau_knn[,c(1,9)]
cau_knn$gene_id <- as.character(cau_knn$gene_id)

ach_arsk_hrsk_cr <- inner_join(ach_arsk_hrsk, cau_rf, by = c('gene' = 'gene_id'))
ach_arsk_hrsk_crs <- inner_join(ach_arsk_hrsk_cr, cau_svr, by = c('gene' = 'gene_id'))
ach_arsk_hrsk_crsk <- inner_join(ach_arsk_hrsk_crs, cau_knn, by = c('gene' = 'gene_id'))

mean(ach_arsk_hrsk_crsk$spearman.x)
mean(ach_arsk_hrsk_crsk$spearman.y)
mean(ach_arsk_hrsk_crsk$spearman)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..x)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..y)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..x.x)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..y.y)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..x.x.x)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..y.y.y)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..x.x.x.x)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d..y.y.y.y)
mean(ach_arsk_hrsk_crsk$spearman_yobs_vs_ypred..d.)



svr <- subset(afa_svr_cv, afa_svr_cv$CV_R2 >= 0.1)
svr$Gene_ID <- as.character(svr$Gene_ID)

rf <- subset(afa_rf_cv, afa_rf_cv$CV_R2 >= 0.1)
rf$Gene_ID <- as.character(rf$Gene_ID)

knn <- subset(afa_knn_cv, afa_knn_cv$CV_R2 >= 0.1)
knn$Gene_ID <- as.character(knn$Gene_ID)

en <- subset(afa_en_cv, afa_en_cv$cv_R2_avg >= 0.1)
en$gene_id <- as.character(en$gene_id)

svr_gen <- inner_join(gencode, svr, by = c("gene_id" = "Gene_ID"))
rf_gen <- inner_join(gencode, rf, by = c("gene_id" = "Gene_ID"))
knn_gen <- inner_join(gencode, knn, by = c("gene_id" = "Gene_ID"))

his_en <- read.table(file = "Z:/data/mesa_models/MESA.HIS.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header = TRUE)
his_en <- subset(mesa_his, cv_R2_avg > 0.1)

cau_en <- read.table(file = "Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header = TRUE)
cau_en <- subset(mesa_cau, cv_R2_avg > 0.1)

#elnet
afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
afa_2_mets$gene <- as.character(afa_2_mets$gene)
cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
cau_2_mets$gene <- as.character(cau_2_mets$gene)
his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
his_2_mets$gene <- as.character(his_2_mets$gene)

ac <- inner_join(afa_2_mets, cau_2_mets, by = c('gene' = 'gene'))
ach <- inner_join(ac, his_2_mets, by = c('gene' = 'gene'))
mean(ach$spearman.x)
mean(ach$spearman.y)
mean(ach$spearman)


#HIS and CAU CV
h_knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/HIS_knn_cv_chr1.txt", header = T, sep = "\t")

#join all chromosomes to create one file for model results of each population
"%&%" <- function(a,b) paste(a,b, sep='')
pop <- "CAU"
mod <- "rf"
model_summaries <- NULL
for (chrom in 1:22) {
  model_summaries <- rbind(model_summaries, read.table(file = "Z:data/mesa_models/python_ml_models/results/" %&% pop %&% "_" %&% mod %&% "_cv_chr" %&% as.character(chrom) %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#Just write out the model_summaries
write.table(model_summaries, file = "Z:/data/mesa_models/python_ml_models/results/" %&% pop %&% "_" %&% mod %&% "_cv_full_chr.txt", row.names = FALSE, quote = F)


# HIS CV
h_knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/HIS_knn_cv_full_chr.txt", header = T)
h_knn <- h_knn[,c(1,2,3)]
h_knn$Gene_Name <- as.character(h_knn$Gene_Name)
h_knn <- subset(h_knn, h_knn$CV_R2 >= 0.1)

h_svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/HIS_svr_rbf_cv_full_chr.txt", header = T)
h_svr <- h_svr[,c(1,2,3)]
h_svr$Gene_Name <- as.character(h_svr$Gene_Name)
h_svr <- subset(h_svr, h_svr$CV_R2 >= 0.1)

h_rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/HIS_rf_cv_full_chr.txt", header = T)
h_rf <- h_rf[,c(1,2,3)]
h_rf$Gene_Name <- as.character(h_rf$Gene_Name)
h_rf <- subset(h_rf, h_rf$CV_R2 >= 0.1)

his_en <- read.table(file = "Z:/data/mesa_models/MESA.HIS.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header = TRUE)
his_en <- his_en[,c(1,2,10)]
his_en$gene_name <- as.character(his_en$gene_name)
his_en <- subset(his_en, his_en$cv_R2_avg >= 0.1)


# CAU CV
c_knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/CAU_knn_cv_full_chr.txt", header = T)

c_svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/CAU_svr_rbf_cv_full_chr.txt", header = T)
c_svr <- c_svr[,c(1,2,3)]
c_svr$Gene_Name <- as.character(c_svr$Gene_Name)
c_svr <- subset(c_svr, c_svr$CV_R2 >= 0.1)

c_rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/CAU_rf_cv_full_chr.txt", header = T)
c_rf <- c_rf[,c(1,2,3)]
c_rf$Gene_Name <- as.character(c_rf$Gene_Name)
c_rf <- subset(c_rf, c_rf$CV_R2 >= 0.1)

cau_en <- read.table(file = "Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header = TRUE)
cau_en <- cau_en[,c(1,2,10)]
cau_en$gene_name <- as.character(cau_en$gene_name)
cau_en <- subset(cau_en, cau_en$cv_R2_avg >= 0.7)
#Check the genes with top performance for each population and model. That is genes at the top 10 or genes with R2 > 0.7

cg <- inner_join(afa_en_cv, afa_rf_cv, by = c("gene_name" = "gene_name")) #AFA:EN+RF
cg = inner_join(cg, afa_svr_cv, by = c("gene_name" = "gene_name")) #AFA:EN+RF+SVR
cg = inner_join(cg, his_en, by = c("gene_name" = "gene_name")) #AFA:EN+RF+SVR+HIS:EN
cg = inner_join(cg, h_rf, by = c("gene_name" = "Gene_Name")) #AFA:EN+RF+SVR+HIS:EN+RF
cg = inner_join(cg, h_svr, by = c("gene_name" = "Gene_Name")) #AFA:EN+RF+SVR+HIS:EN+RF+SVR
cg = inner_join(cg, cau_en, by = c("gene_name" = "gene_name")) #AFA:EN+RF+SVR+HIS:EN+RF+SVR+CAU:EN
cg = inner_join(cg, c_rf, by = c("gene_name" = "Gene_Name")) #AFA:EN+RF+SVR+HIS:EN+RF+SVR+CAU:EN+RF
cg = inner_join(cg, c_svr, by = c("gene_name" = "Gene_Name")) #AFA:EN+RF+SVR+HIS:EN+RF+SVR+CAU:EN+RF+SVR


#MESA 2 METS
#Change the gencode to remove the decimal in gene_id
gnc <- gencode
for (i in 1:length(gnc$gene_id)){
  gnc$gene_id[i] <- gsub('\\.[0-9]+','',gnc$gene_id[i])
}

#elastic net
a2m <- inner_join(gnc, afa_2_mets, by = c("gene_id" = "gene"))
c2m <- inner_join(gnc, cau_2_mets, by = c("gene_id" = "gene"))
h2m <- inner_join(gnc, his_2_mets, by = c("gene_id" = "gene"))

#Random Forest
a2m_rf <- afa_rf
a2m_rf <- inner_join(gnc, a2m_rf, by = c("gene_id" = "gene_id"))

c2m_rf <- cau_rf
c2m_rf <- inner_join(gnc, c2m_rf, by = c("gene_id" = "gene_id"))

h2m_rf <- his_rf
h2m_rf <- inner_join(gnc, h2m_rf, by = c("gene_id" = "gene_id"))

#SVM
a2m_svr <- afa_svr
a2m_svr <- inner_join(gnc, a2m_svr, by = c("gene_id" = "gene_id"))

c2m_svr <- cau_svr
c2m_svr <- inner_join(gnc, c2m_svr, by = c("gene_id" = "gene_id"))

h2m_svr <- his_svr
h2m_svr <- inner_join(gnc, h2m_svr, by = c("gene_id" = "gene_id"))


#For approximately what percent of genes with rho>0.1, is RF rho higher than EN?
#FOR AFA
afa2mets <- inner_join(afa_2_mets, afa_rf, by = c("gene" = "gene_id"))
names(afa2mets) <- c("genes", "EN", "RF")
#afa2mets <- subset(afa2mets, EN > 0.1)# & RF > 0.1)

en <- 0 #Count how many elastic net genes have rho > 0.1 & > rf genes
rf <- 0 #Count how many random forest genes have rho > 0.1 & > EN genes

for (i in 1:length(afa2mets$genes)){
  if (afa2mets$EN[i] > 0.1 & afa2mets$EN[i] > afa2mets$RF[i]){
    en <- en + 1
  }
  if (afa2mets$RF[i] > 0.1 & afa2mets$RF[i] > afa2mets$EN[i]){
    rf <- rf + 1
  }
}

#FOR HIS
his2mets <- inner_join(his_2_mets, his_rf, by = c("gene" = "gene_id"))
names(his2mets) <- c("genes", "EN", "RF")
#afa2mets <- subset(afa2mets, EN > 0.1)# & RF > 0.1)

en <- 0 #Count how many elastic net genes have rho > 0.1 & > rf genes
rf <- 0 #Count how many random forest genes have rho > 0.1 & > EN genes

for (i in 1:length(his2mets$genes)){
  if (his2mets$EN[i] > 0.1 & his2mets$EN[i] > his2mets$RF[i]){
    en <- en + 1
  }
  if (his2mets$RF[i] > 0.1 & his2mets$RF[i] > his2mets$EN[i]){
    rf <- rf + 1
  }
}


#FOR CAU
cau2mets <- inner_join(cau_2_mets, cau_rf, by = c("gene" = "gene_id"))
names(cau2mets) <- c("genes", "EN", "RF")
#afa2mets <- subset(afa2mets, EN > 0.1)# & RF > 0.1)

en <- 0 #Count how many elastic net genes have rho > 0.1 & > rf genes
rf <- 0 #Count how many random forest genes have rho > 0.1 & > EN genes

for (i in 1:length(cau2mets$genes)){
  if (cau2mets$EN[i] > 0.1 & cau2mets$EN[i] > cau2mets$RF[i]){
    en <- en + 1
  }
  if (cau2mets$RF[i] > 0.1 & cau2mets$RF[i] > cau2mets$EN[i]){
    rf <- rf + 1
  }
}
