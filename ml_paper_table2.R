#Elastic Net
afhi_2_mets_en <- read.table(file="Z:/data/ml_paper/afhi_2_mets_en_rho_all.txt", header=T)
all_2_mets_en <- read.table(file="Z:/data/ml_paper/all_2_mets_en_rho_all.txt", header=T)
cau_2_mets_en <- read.table(file="Z:/data/ml_paper/cau_2_mets_en_rho_all.txt", header=T)
his_2_mets_en <- read.table(file="Z:/data/ml_paper/his_2_mets_en_rho_all.txt", header=T)
afa_2_mets_en <- read.table(file="Z:/data/ml_paper/afa_2_mets_en_rho_all.txt", header=T)

#take only genes in common across population
library(dplyr)
en_merge_pop <- inner_join(afhi_2_mets_en, all_2_mets_en, by = c("gene"="gene"))
en_merge_pop <- inner_join(en_merge_pop, cau_2_mets_en, by = c("gene"="gene"))
en_merge_pop <- inner_join(en_merge_pop, his_2_mets_en, by = c("gene"="gene"))
en_merge_pop <- inner_join(en_merge_pop, afa_2_mets_en, by = c("gene"="gene"))
names(en_merge_pop) <- c("gene", "afhi", "all", "cau", "his", "afa")






#Random Forest
afhi_2_mets_rf <- read.table(file="Z:/data/ml_paper/afhi_2_mets_rf_rho_all.txt", header=T)
all_2_mets_rf <- read.table(file="Z:/data/ml_paper/all_2_mets_rf_rho_all.txt", header=T)
cau_2_mets_rf <- read.table(file="Z:/data/ml_paper/cau_2_mets_rf_rho_all.txt", header=T)
his_2_mets_rf <- read.table(file="Z:/data/ml_paper/his_2_mets_rf_rho_all.txt", header=T)
afa_2_mets_rf <- read.table(file="Z:/data/ml_paper/afa_2_mets_rf_rho_all.txt", header=T)

#take only grfes in common across population
library(dplyr)
rf_merge_pop <- inner_join(afhi_2_mets_rf, all_2_mets_rf, by = c("gene"="gene"))
rf_merge_pop <- inner_join(rf_merge_pop, cau_2_mets_rf, by = c("gene"="gene"))
rf_merge_pop <- inner_join(rf_merge_pop, his_2_mets_rf, by = c("gene"="gene"))
rf_merge_pop <- inner_join(rf_merge_pop, afa_2_mets_rf, by = c("gene"="gene"))
names(rf_merge_pop) <- c("gene", "afhi", "all", "cau", "his", "afa")




#Support Vector
afhi_2_mets_svr <- read.table(file="Z:/data/ml_paper/afhi_2_mets_svr_rho_all.txt", header=T)
all_2_mets_svr <- read.table(file="Z:/data/ml_paper/all_2_mets_svr_rho_all.txt", header=T)
cau_2_mets_svr <- read.table(file="Z:/data/ml_paper/cau_2_mets_svr_rho_all.txt", header=T)
his_2_mets_svr <- read.table(file="Z:/data/ml_paper/his_2_mets_svr_rho_all.txt", header=T)
afa_2_mets_svr <- read.table(file="Z:/data/ml_paper/afa_2_mets_svr_rho_all.txt", header=T)

#take only gsvres in common across population
library(dplyr)
svr_merge_pop <- inner_join(afhi_2_mets_svr, all_2_mets_svr, by = c("gene"="gene"))
svr_merge_pop <- inner_join(svr_merge_pop, cau_2_mets_svr, by = c("gene"="gene"))
svr_merge_pop <- inner_join(svr_merge_pop, his_2_mets_svr, by = c("gene"="gene"))
svr_merge_pop <- inner_join(svr_merge_pop, afa_2_mets_svr, by = c("gene"="gene"))
names(svr_merge_pop) <- c("gene", "afhi", "all", "cau", "his", "afa")





# K Nearest Neighbour
afhi_2_mets_knn <- read.table(file="Z:/data/ml_paper/afhi_2_mets_knn_rho_all.txt", header=T)
all_2_mets_knn <- read.table(file="Z:/data/ml_paper/all_2_mets_knn_rho_all.txt", header=T)
cau_2_mets_knn <- read.table(file="Z:/data/ml_paper/cau_2_mets_knn_rho_all.txt", header=T)
his_2_mets_knn <- read.table(file="Z:/data/ml_paper/his_2_mets_knn_rho_all.txt", header=T)
afa_2_mets_knn <- read.table(file="Z:/data/ml_paper/afa_2_mets_knn_rho_all.txt", header=T)

#take only gknnes in common across population
library(dplyr)
knn_merge_pop <- inner_join(afhi_2_mets_knn, all_2_mets_knn, by = c("gene"="gene"))
knn_merge_pop <- inner_join(knn_merge_pop, cau_2_mets_knn, by = c("gene"="gene"))
knn_merge_pop <- inner_join(knn_merge_pop, his_2_mets_knn, by = c("gene"="gene"))
knn_merge_pop <- inner_join(knn_merge_pop, afa_2_mets_knn, by = c("gene"="gene"))
names(knn_merge_pop) <- c("gene", "afhi", "all", "cau", "his", "afa")
