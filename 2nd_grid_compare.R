c11 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_rf_cv_chr1_chunk1.txt", header = T)
c12 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_rf_cv_chr1_chunk2.txt", header = T)
c13 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_rf_cv_chr1_chunk3.txt", header = T)
c14 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_rf_cv_chr1_chunk4.txt", header = T)
c15 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_rf_cv_chr1_chunk5.txt", header = T)



an <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_chr1_chunk3.txt", header = T)
ap <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_parameter_per_gene_chr1_chunk3.txt", header = T)
chr <- read.table(file = "Z:/data/mesa_models/python_ml_models/laptop_results/2nd_best_grid_split_rf_cv_chr6_chunk2.txt", header = T)


oc11 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_rf_cv_chr1_chunk1.txt", header = T)
oc12 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_rf_cv_chr1_chunk2.txt", header = T)
oc13 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_rf_cv_chr1_chunk3.txt", header = T)
oc14 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_rf_cv_chr1_chunk4.txt", header = T)
oc15 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_rf_cv_chr1_chunk5.txt", header = T)


oc12knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_knn_cv_chr1_chunk2.txt", header = T)
c12knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_knn_cv_chr1_chunk2.txt", header = T)

oc12svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_svr_cv_chr1_chunk2.txt", header = T)
c12ksvr <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_svr_cv_chr1_chunk2.txt", header = T)

oc12el <- read.table(file = "Z:/data/mesa_models/split_mesa/results/full_AFA_chr1_model_summaries.txt", header = T)

# merge all the chunk results
c62 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_chr6_chunk2.txt", header = T)
c62 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_parameter_per_gene_chr6_chunk2.txt", header = T)
c62 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_svr_grid_split_parameter_per_gene_chr6_chunk2.txt", header = T, sep = "\t")
c62 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_knn_grid_split_parameter_per_gene_chr6_chunk2.txt", header = T)

#step 1
"%&%" <- function(a,b) paste(a,b, sep = "")


#alg <- "knn"

#Merging only the old chunks that were complete
algs <- c("rf", "svr", "knn")
no <- 21

for (alg in algs){
  oc21 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk1.txt", header = T)
  oc22 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk2.txt", header = T)
  oc23 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk3.txt", header = T)
  oc24 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk4.txt", header = T)
  oc25 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk5.txt", header = T)
  
  full <- rbind(oc21, oc22, oc23, oc24, oc25)
  write.table(full, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_chr"%&% no %&% "_full.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  
}


#Merging the old and new chunks that are complete
"%&%" <- function(a,b) paste(a,b, sep = "")

alg <- "knn"
#no <- 1

algs <- c("rf", "svr", "knn")
no <- 22

for (alg in algs){
  oc21 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk1.txt", header = T)
  oc22 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk2.txt", header = T)
  oc23 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk3.txt", header = T)
  oc24 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk4.txt", header = T)
  oc25 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk5.txt", header = T)
  
  c21 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk1.txt", header = T)
  c22 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk2.txt", header = T)
  c23 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk3.txt", header = T)
  c24 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk4.txt", header = T)
  c25 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk5.txt", header = T)
  
  full <- rbind(oc21, oc22, oc23, oc24, oc25, c21, c22, c23, c24, c25)
  #full <- rbind(oc21, oc22, oc23, oc24, oc25, c22, c24, c25)
  full <- full[unique(full$Gene_ID),]
  write.table(full, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_chr"%&% no %&% "_full.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  
}

octest <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_chr"%&% no %&% "_full.txt", header = T)
alg <- 'rf'
no <- 22
oc21 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk1.txt", header = T)
oc22 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk2.txt", header = T)
oc23 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk3.txt", header = T)
oc24 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk4.txt", header = T)
oc25 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk5.txt", header = T, sep = "\t")

c21 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk1.txt", header = T)
c22 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk2.txt", header = T)
c23 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk3.txt", header = T)
c24 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk4.txt", header = T)
c25 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk5.txt", header = T)

full <- rbind(oc21, oc22, oc23, oc24, oc25, c21, c22, c23, c24, c25)
full <- full[unique(full$Gene_ID),]
write.table(full, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_chr"%&% no %&% "_full.txt", quote = FALSE, sep = "\t", row.names = FALSE)

octest <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_chr"%&% no %&% "_full.txt", header = T)


#Step 2
# merge all chromosome grid to one
"%&%" <- function(a,b) paste(a,b, sep = "")
algs <- c("rf", "svr", "knn")

comp <- NULL

for (alg in algs){
  for (no in 1:22){
    comp <- rbind(comp, read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&%"_chr"%&% no %&%"_full.txt", header = T))
  }
  write.table(comp, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_all_chrom.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  
}

#write.table(comp, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_full_chrom.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#merge all chromosomes manually
"%&%" <- function(a,b) paste(a,b, sep = "")
rf_comp <- NULL
alg <- "rf"
for (no in 1:22){
  rf_comp <- rbind(rf_comp, read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&%"_chr"%&% no %&%"_full.txt", header = T))
}
write.table(rf_comp, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_all_chrom.txt", quote = FALSE, sep = "\t", row.names = FALSE)


svr_comp <- NULL
alg <- "svr"
for (no in 1:22){
  svr_comp <- rbind(svr_comp, read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&%"_chr"%&% no %&%"_full.txt", header = T))
}
write.table(svr_comp, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_all_chrom.txt", quote = FALSE, sep = "\t", row.names = FALSE)


knn_comp <- NULL
alg <- "knn"
for (no in 1:22){
  knn_comp <- rbind(knn_comp, read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&%"_chr"%&% no %&%"_full.txt", header = T))
}
write.table(knn_comp, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_"%&% alg %&% "_all_chrom.txt", quote = FALSE, sep = "\t", row.names = FALSE)





# Chapter 2
# compare the optimized machine learning models (rf, svr, knn) with R elastic net in predixcan model

#full elnet AFA chrom
elnet <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
rf_com <- read.table(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
#Random forest

#sort by gene_id
rf_comp$Gene_ID <- as.character(rf_comp$Gene_ID)
elnet$gene_id <- as.character(elnet$gene_id)

#sort the dataframe by gene id
rf_comp <- rf_comp[order(rf_comp$Gene_ID),]
elnet <- elnet[order(elnet$gene_id),]

#before filtering, they both have 9623 genes

#filter the R2 to be within -1 and 1
rf_filt <- subset(rf_comp, CV_R2 >= -1)
elnet_filt <- subset(elnet, cv_R2_avg >= -1)
# after filtering, random forest has 9622, while elnet has 9609 genes left

rf_filt <- rf_filt[,c(1,3)]
elnet_filt <- elnet_filt[,c(1,10)]

rf_filt$Gene_ID <- as.character(rf_filt$Gene_ID)
elnet_filt$gene_id <- as.character(elnet_filt$gene_id)

library(dplyr)

el_rf_filt <- inner_join(elnet_filt, rf_filt, by = c("gene_id" = "Gene_ID"))
names(el_rf_filt) <- c("gene_id", "elnet_cv_R2", "rf_cv_R2")

#hist(el_rf_filt, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Baseline (Elnet) Performance on MESA AFA all chr")
#hist(rf_comp$CV_R2, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", main = "Random Forest Performance on MESA AFA all chr")

ggplot(el_rf_filt, aes(x=elnet_cv_R2, y=rf_cv_R2)) + 
  ggtitle("Elastic Net vs Random Forest (outliers (< -1) filtered off) Grid Search") + 
  ylab("RF R2 Grid Optimized") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

#check the correlation of the two: elnet_filt and rf_filt

library("ggpubr")
ggscatter(el_rf_filt, x = "elnet_cv_R2", y = "rf_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net R2", ylab = "RF R2 Grid Optimized", 
          title = "Elastic Net vs Random Forest (outliers (< -1) filtered off) Grid Search",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR
#sort the dataframe by gene id
svr_comp <- svr_comp[order(svr_comp$Gene_ID),]
elnet <- elnet[order(elnet$gene_id),]

#filter the R2 to be within -1 and 1

svr_filt <- subset(svr_comp, CV_R2 >= -1)


svr_filt <- svr_filt[,c(1,3)]
elnet_filt <- elnet_filt[,c(1,10)]

svr_filt$Gene_ID <- as.character(svr_filt$Gene_ID)
elnet_filt$gene_id <- as.character(elnet_filt$gene_id)

library(dplyr)

el_svr_filt <- inner_join(elnet_filt, svr_filt, by = c("gene_id" = "Gene_ID"))
names(el_svr_filt) <- c("gene_id", "elnet_cv_R2", "svr_cv_R2")

ggplot(el_svr_filt, aes(x=elnet_cv_R2, y=svr_cv_R2)) + 
  ggtitle("Elastic Net vs SVM (outliers (< -1) filtered off) Grid Search") + 
  ylab("SVM R2 Grid Optimized") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

#check the correlation of the two: elnet_filt and rf_filt

library("ggpubr")
ggscatter(el_svr_filt, x = "elnet_cv_R2", y = "svr_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net R2", ylab = "SVM R2 Grid Optimized", 
          title = "Elastic Net vs SVM (outliers (< -1) filtered off) Grid Search",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#KNN
#sort the dataframe by gene id
knn_comp <- knn_comp[order(knn_comp$Gene_ID),]
elnet <- elnet[order(elnet$gene_id),]

#filter the R2 to be within -1 and 1

knn_filt <- subset(knn_comp, CV_R2 >= -1)


knn_filt <- knn_filt[,c(1,3)]
elnet_filt <- elnet_filt[,c(1,10)]

knn_filt$Gene_ID <- as.character(knn_filt$Gene_ID)
elnet_filt$gene_id <- as.character(elnet_filt$gene_id)

library(dplyr)

el_knn_filt <- inner_join(elnet_filt, knn_filt, by = c("gene_id" = "Gene_ID"))
names(el_knn_filt) <- c("gene_id", "elnet_cv_R2", "knn_cv_R2")

ggplot(el_knn_filt, aes(x=elnet_cv_R2, y=knn_cv_R2)) + 
  ggtitle("Elastic Net vs KNN (outliers (< -1) filtered off) Grid Search") + 
  ylab("KNN R2 Grid Optimized") + xlab("Elastic Net R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

#check the correlation of the two: elnet_filt and rf_filt

library("ggpubr")
ggscatter(el_knn_filt, x = "elnet_cv_R2", y = "knn_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net R2", ylab = "KNN R2 Grid Optimized", 
          title = "Elastic Net vs SVM (outliers (< -1) filtered off) Grid Search",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#Chapter 3
#compare the grid optimized model with the un-optimized ones
#read in saved files
#RF
rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/rf_cv_full_AFA_chr.txt", header = TRUE)

rf_comp <- rf_comp[order(rf_comp$Gene_ID),]
rf <- rf[order(rf$Gene_ID),]

#filter the R2 to be within -1 and 1

knn_filt <- subset(knn_comp, CV_R2 >= -1)


orf <- rf_comp[,c(1,3)]
rf <- rf[,c(1,2)]

orf$Gene_ID <- as.character(orf$Gene_ID)
rf$Gene_ID <- as.character(rf$Gene_ID)

library(dplyr)

orf_rf <- inner_join(orf, rf, by = c("Gene_ID" = "Gene_ID"))
names(orf_rf) <- c("gene_id", "orf_cv_R2", "rf_cv_R2")

ggplot(orf_rf, aes(x=orf_cv_R2, y=rf_cv_R2)) + 
  ggtitle("optimized RF vs non optimized RF") + 
  ylab("Non Optimized RF") + xlab("Optimized RF R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

#check the correlation of the two: elnet_filt and rf_filt

library("ggpubr")
ggscatter(orf_rf, x = "orf_cv_R2", y = "rf_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized RF R2", ylab = "Non Optimized RF", 
          title = "optimized RF vs non optimized RF",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#RF
knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/knn_cv_full_AFA_chr.txt", header = TRUE)

knn_comp <- knn_comp[order(knn_comp$Gene_ID),]
knn <- knn[order(knn$Gene_ID),]

#filter the R2 to be within -1 and 1

#knn_filt <- subset(knn_comp, CV_R2 >= -1)


oknn <- knn_comp[,c(1,3)]
knn <- knn[,c(1,2)]

oknn$Gene_ID <- as.character(oknn$Gene_ID)
knn$Gene_ID <- as.character(knn$Gene_ID)

library(dplyr)

oknn_knn <- inner_join(oknn, knn, by = c("Gene_ID" = "Gene_ID"))
names(oknn_knn) <- c("gene_id", "oknn_cv_R2", "knn_cv_R2")

ggplot(oknn_knn, aes(x=oknn_cv_R2, y=knn_cv_R2)) + 
  ggtitle("optimized KNN vs non optimized KNN") + 
  ylab("Non Optimized KNN") + xlab("Optimized KNN R2") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1, 1)) + ylim(c(-1,1))

#check the correlation of the two: elnet_filt and rf_filt

library("ggpubr")
ggscatter(oknn_knn, x = "oknn_cv_R2", y = "knn_cv_R2", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized KNN R2", ylab = "Non Optimized KNN", 
          title = "optimized KNN vs non optimized KNN",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


mets1 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_rf_cor_test_chr22.txt", header = T, sep = "\t")
omets1 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_rf_cor_test_chr22.txt", header = T, sep = "\t")


gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
gv28 <- read.table(file = "Z:/data/METS_model/hg19/gencode.v28_annotation.parsed.txt", header = T)


rf_grid <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_chr20_full.txt", header=T, sep="\t")

afa_g_no <- read.table(file = "Z:/data/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
mets_g_no <- read.table(file = "Z:/data/METS_model/hg19/METS_peer10_all_chr_prot_coding_gex.txt", header = T)

library(tidyverse)

afa_g_no$PROBE_ID <- as.character(afa_g_no$PROBE_ID)
mets_g_no$PROBE_ID <- as.character(mets_g_no$PROBE_ID)

no = 1
chrom <- subset(gv18, gv18$chr == no)
chrom <- subset(chrom, chrom$gene_type == "protein_coding")
chrom$gene_id <- as.character(chrom$gene_id)
g_ex <- inner_join(chrom, afa_g_no, by = c("gene_id" = "PROBE_ID"))

mchrom <- subset(gv28, gv28$chr == no)
mchrom <- subset(mchrom, mchrom$gene_type == "protein_coding")
mchrom$gene_id <- as.character(mchrom$gene_id)
mg_ex <- inner_join(mchrom, mets_g_no, by = c("gene_id" = "PROBE_ID"))

afa <- g_ex
met <- mg_ex

#am <- inner_join(afa[,2], met[,2], by = c("gene_id" = "gene_id"))

for (i in 1:length(afa$gene_id)){
  afa$gene_id[i] <- gsub('\\.[0-9]+','',afa$gene_id[i])
} #just to remove the decimal places in the gene_id


for (i in 1:length(met$gene_id)){
  met$gene_id[i] <- gsub('\\.[0-9]+','',met$gene_id[i])
} #just to remove the decimal places in the gene_id

uni <- afa$gene_id %in% met$gene_id
uni <- uni[uni != FALSE]

a2m9 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_rf_cor_test_chr20.txt", header = T, sep = "\t")
oa2m9 <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_rf_cor_test_chr20.txt", header = T, sep = "\t")


# Model performance on test data set

"%&%" <- function(a,b) paste(a,b, sep = "")

knn <- NULL
rf <- NULL
svr <- NULL
#el <- NULL
pop <- "AFA"

for (chrom in 1:22) {
  no <- as.character(chrom)
  knn <- rbind(knn, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_knn_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  rf <- rbind(rf, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_rf_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  svr <- rbind(svr, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_svr_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(knn, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_knn_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rf, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_rf_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_svr_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#Compare model Optimized performance on AFA 2 METS as against Elastic Net
#AFA 2METS

afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
afa_2_mets$gene <- as.character(afa_2_mets$gene)

knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_knn_cor_test_full_chr.txt", header = T)
knn <- knn[,c(1,9)]
knn$gene_id <- as.character(knn$gene_id)

rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_rf_cor_test_full_chr.txt", header = T)
rf <- rf[,c(1,9)]
rf$gene_id <- as.character(rf$gene_id)

svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_svr_cor_test_full_chr.txt", header = T)
svr <- svr[,c(1,9)]
svr$gene_id <- as.character(svr$gene_id)

#Take the overlapping genes between elastic net and each model
#RF
library(dplyr)
library(ggplot2)
library("ggpubr")
elrf <- inner_join(afa_2_mets, rf, by = c("gene" = "gene_id"))
elrf <- elrf[,c(2,3)]
names(elrf) <- c("elnet", "rf")

#Before plotting, take the difference of rf and elnet (rf - elnet) plot on y axis against elnet on x axis
rfdif <- elrf$rf - elrf$elnet

ggplot(elrf, aes(x=elnet, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("RF Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elrf, x = "elnet", y = "rf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "RF Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs (rf - elnet)
elrf <- cbind(elrf, rfdif)

ggplot(elrf, aes(x=elnet, y=rfdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("RF - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where knn minus elnet > 0. That is, how many genes where svr outperformed elnet. total genes = 3227
sum(elrf$rfdif > 0) # answer = 1443.Therefore elnet performed better only on 1784 genes

#Compare with non optimized and optimized rf
orf <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_rf_cor_test_all_chr.txt", header = T, sep = "\t")
orf <- orf[,c(1,9)]
orf$gene_id <- as.character(orf$gene_id)

rforf <- inner_join(rf, orf, by = c("gene_id" = "gene_id"))
names(rforf) <- c("gene", "rf", "orf")

ggplot(rforf, aes(x=rf, y=orf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("Non Optimized RF") + xlab("Optimized RF") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(rforf, x = "rf", y = "orf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized RF", ylab = "Non Optimized RF", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#SVR
library(dplyr)
library(ggplot2)
library("ggpubr")
elsvr <- inner_join(afa_2_mets, svr, by = c("gene" = "gene_id"))
elsvr <- elsvr[,c(2,3)]
names(elsvr) <- c("elnet", "svr")
ggplot(elsvr, aes(x=elnet, y=svr)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("SVR Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elsvr, x = "elnet", y = "svr", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs SVR - Elnet
svrdif <- elsvr$svr - elsvr$elnet
elsvr <- cbind(elsvr, svrdif)

ggplot(elsvr, aes(x=elnet, y=svrdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("SVR - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where svr minus elnet > 0. That is, how many genes where svr outperformed elnet. total genes = 3227
sum(elsvr$svrdif > 0) # answer = 1297.Therefore elnet performed better only on 1930 genes

#Compare with non optimized and optimized svr
osvr <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_svr_rbf_cor_test_all_chr.txt", header = T, sep = "\t")
osvr <- osvr[,c(1,9)]
osvr$gene_id <- as.character(osvr$gene_id)

#KNN
library(dplyr)
library(ggplot2)
library("ggpubr")
elknn <- inner_join(afa_2_mets, knn, by = c("gene" = "gene_id"))
elknn <- elknn[,c(2,3)]
names(elknn) <- c("elnet", "knn")
ggplot(elknn, aes(x=elnet, y=knn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("KNN Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elknn, x = "elnet", y = "knn", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs knn - elnet
knndif <- elknn$knn - elknn$elnet
elknn <- cbind(elknn, knndif)
ggplot(elknn, aes(x=elnet, y=knndif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("KNN - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where knn minus elnet > 0. That is, how many genes where svr outperformed elnet. total genes = 3228
sum(elknn$knndif > 0) # answer = 1221.Therefore elnet performed better only on 2007 genes

#Compare with non optimized and optimized KNN
oknn <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_knn_cor_test_all_chr.txt", header = T, sep = "\t")
oknn <- oknn[,c(1,9)]
oknn$gene_id <- as.character(oknn$gene_id)

knnoknn <- inner_join(knn, oknn, by = c("gene_id" = "gene_id"))
names(knnoknn) <- c("gene", "knn", "oknn")

ggplot(knnoknn, aes(x=knn, y=oknn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("Non Optimized KNN") + xlab("Optimized KNN") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(knnoknn, x = "knn", y = "oknn", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized KNN", ylab = "Non Optimized KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


###############################
#Random Forest Trees
trees <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
c1c1t <- read.table(file="Z:/data/mesa_models/python_ml_models/results/grid_split/AFA_rf_grid_split_chr1_chunk1.txt", header=T)

#Just plot first gene 
g1 <- as.numeric(c1c1t[1,c(4:13)])
#g1 <- g1[1]

plot(trees, g1)


#############################
#KNN
k <- c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)
c1c1k <- read.table(file="Z:/data/mesa_models/python_ml_models/results/grid_split/AFA_knn_grid_split_parameter_per_gene_chr1_chunk1.txt", header=T, sep="\t")



#Check the performance of Random Forest on 5 trees only
#AFA CV R2
#Merge Chromosomes

rf5 <- read.table(file="Z:/data/mesa_models/python_ml_models/results/CAU_rf_tree_check_cv_chr1.txt", header=T, sep="\t")


"%&%" <- function(a,b) paste(a,b, sep = "")

rf_afa <- NULL
pop <- "AFA"

for (chrom in 1:22) {
  no <- as.character(chrom)
  rf_afa <- rbind(rf_afa, read.table(file = "Z:/data/mesa_models/python_ml_models/results/" %&% pop %&% "_rf_tree_check_cv_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(rf_afa, file = "Z:/data/mesa_models/python_ml_models/results/" %&% pop %&% "_rf_5tree_cv_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#HIS CV R2
"%&%" <- function(a,b) paste(a,b, sep = "")

rf_his <- NULL
pop <- "HIS"

for (chrom in 1:22) {
  no <- as.character(chrom)
  rf_his <- rbind(rf_his, read.table(file = "Z:/data/mesa_models/python_ml_models/results/" %&% pop %&% "_rf_tree_check_cv_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(rf_his, file = "Z:/data/mesa_models/python_ml_models/results/" %&% pop %&% "_rf_5tree_cv_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#CAU CV R2
"%&%" <- function(a,b) paste(a,b, sep = "")

rf_cau <- NULL
pop <- "CAU"

for (chrom in 1:22) {
  no <- as.character(chrom)
  rf_cau <- rbind(rf_cau, read.table(file = "Z:/data/mesa_models/python_ml_models/results/" %&% pop %&% "_rf_tree_check_cv_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(rf_cau, file = "Z:/data/mesa_models/python_ml_models/results/" %&% pop %&% "_rf_5tree_cv_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#compare the optimized RF with the 5 tree
#AFA
#elnet <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
rf_afa_op <- read.table(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
rf_afa_op <- rf_afa_op[,c(1,3)]
rf_afa_op$Gene_ID <- as.character(rf_afa_op$Gene_ID)
rf_afa_op <- subset(rf_afa_op, CV_R2 >= -1)

rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/AFA_rf_5tree_cv_full_chr.txt", header=T, sep="\t")
rf_afa <- rf_afa[,c(1,3)]
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- subset(rf_afa, CV_R2 >= -1)

library(dplyr)
library(ggpubr)

op_5tr_afa <- inner_join(rf_afa_op, rf_afa, by = c("Gene_ID" = "Gene_ID"))
names(op_5tr_afa) <- c('gene', 'op', 'tr5')

ggplot(op_5tr_afa, aes(x=op, y=tr5)) + ggtitle("AFA CV R2 Optimized RF vs 5 Trees RF") + 
  ylab("5 Trees RF") + xlab("Optimized RF") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(op_5tr_afa, x = "op", y = "tr5", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized RF", ylab = "5 Trees RF", title = "AFA CV R2 Optimized RF vs 5 Trees RF",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#HIS
rf_afa_op <- read.table(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
rf_afa_op <- rf_afa_op[,c(1,3)]
rf_afa_op$Gene_ID <- as.character(rf_afa_op$Gene_ID)
rf_afa_op <- subset(rf_afa_op, CV_R2 >= -1)

rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/AFA_rf_5tree_cv_full_chr.txt", header=T, sep="\t")
rf_afa <- rf_afa[,c(1,3)]
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf_afa <- subset(rf_afa, CV_R2 >= -1)

library(dplyr)
library(ggpubr)

op_5tr_afa <- inner_join(rf_afa_op, rf_afa, by = c("Gene_ID" = "Gene_ID"))
names(op_5tr_afa) <- c('gene', 'op', 'tr5')

ggplot(op_5tr_afa, aes(x=op, y=tr5)) + ggtitle("AFA CV R2 Optimized RF vs 5 Trees RF") + 
  ylab("5 Trees RF") + xlab("Optimized RF") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(op_5tr_afa, x = "op", y = "tr5", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized RF", ylab = "5 Trees RF", title = "AFA CV R2 Optimized RF vs 5 Trees RF",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")


#AFA2METS 5trees

"%&%" <- function(a,b) paste(a,b, sep = "")

a2me <- NULL
pop <- "AFA"

for (chrom in 1:22) {
  no <- as.character(chrom)
  a2me <- rbind(a2me, read.table(file = "Z:/data/mesa_models/python_ml_models/results/rf_5tree_" %&% pop %&% "_2_METS_cor_test_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(a2me, file = "Z:/data/mesa_models/python_ml_models/results/rf_5tree_" %&% pop %&% "_2_METS_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)

a2me <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/rf_5tree_AFA_2_METS_cor_test_full_chr.txt", header=T, sep="\t")
a2me <- a2me[,c(1,9)]
a2me$gene_id <- as.character(a2me$gene_id)

rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_rf_cor_test_full_chr.txt", header = T)
rf <- rf[,c(1,9)]
rf$gene_id <- as.character(rf$gene_id)

a2m_op_5tr <- inner_join(rf, a2me, by = c("gene_id" = "gene_id"))
names(a2m_op_5tr) <- c('gene', 'op', 'tr5')

ggplot(a2m_op_5tr, aes(x=op, y=tr5)) + ggtitle("AFA 2 METS cor Optimized RF vs 5 Trees RF") + 
  ylab("5 Trees RF") + xlab("Optimized RF") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(a2m_op_5tr, x = "op", y = "tr5", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized RF", ylab = "5 Trees RF", title = "AFA 2 METS cor Optimized RF vs 5 Trees RF",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")




############################################################################################################
#Merge all the RF tree params
"%&%" <- function(a,b) paste(a,b, sep = "")

prf1 <- NULL

for (chrom in 1:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    prf1 <- rbind(prf1, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_split/AFA_rf_grid_split_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header = T, sep = "\t"))
  }
}

#Merge 2nd Chunk

prf2 <- NULL

for (chrom in c(1,12,17,19,20,22)) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
      prf2 <- rbind(prf2, read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header = T, sep = "\t"))
  }
}

#chrom 6
prf3 <- read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_chr6_chunk2.txt", header = T, sep = "\t")

prf4 <- NULL
no = '14' #chrom 14
for (chunk in 2:5){
    prf4 <- rbind(prf4, read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header = T, sep = "\t"))
}

prf5 <- NULL
no = '16' #chrom 16
for (chunk in c(2,4,5)){
  prf5 <- rbind(prf5, read.table(file = "Z:/data/mesa_models/python_ml_models/new_results/2nd_AFA_rf_grid_split_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header = T, sep = "\t"))
}

#join all
prf <- rbind(prf1,prf2,prf3,prf4,prf5)
prf$X <- NULL
prf$gene_name <- NULL
prf$chr <- NULL

# The Unique function works only when the files are read in without including stringAsFactor=F. Just ignore this

prf <- prf[unique(prf$gene_id),]

write.table(prf, file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_rf_100_to_500tree_cv_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#Plot the different trees of the Random Forest
brf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header=T)
prf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_rf_100_to_500tree_cv_full_chr.txt", header=T)

#Start by taking genes with best CV R2 > 0.3
#Then use the genes to subset the dataframe of random forest trees
brf <- subset(brf, CV_R2 > 0.5)
brf <- brf[,c(1,2)]
brf$Gene_ID <- as.character(brf$Gene_ID)

prf <- prf[,c(1,4:13)]
prf$gene_id <- as.character(prf$gene_id)

library(tidyverse)

rftrees <- inner_join(brf, prf, by = c("Gene_ID" = "gene_id"))
#I had to write out the dataframe to file so as to be able to use the "unique' function to take the unique genes in the dataframe
write.table(rftrees, file = "Z:/data/mesa_models/python_ml_models/rf_cv_gene_param.txt", quote = FALSE, sep = "\t", row.names = FALSE)
rftrees <- read.table(file = "Z:/data/mesa_models/python_ml_models/rf_cv_gene_param.txt", header=T)
#rftrees <- as.data.frame(rftrees)
rftrees <- rftrees[unique(rftrees$Gene_ID),]
names(rftrees) <- c('gene_id', 'gene_name', '50','100','150','200','250','300','350','400','450','500')
rftrees <- rftrees[,c(3:12)]

trees <- c(50,100,150,200,250,300,350,400,450,500)

#matplot(rftrees, type = c("p"),pch=1,col = 1:4) #plot

cl <- rainbow(length(rftrees$`50`))
plot(trees, rftrees[1,], type = 'o', xlab = "Trees", ylab = "CV R2", main = "RF Tree", ylim = c(0.5,1), col=cl[1])

for (i in 2:length(rftrees$`50`)){
  lines(trees, rftrees[i,], type = "o", col=cl[i])
}

library(ggplot2)

p = ggplot(rftrees) + 
  geom_line(aes(x = trees, y = rftrees[1,]), color = "blue")# +
  #geom_line(data = prescription2, aes(x = dates, y = Difference), color = "red") +
  #xlab('Dates') +
  #ylab('percent.change')

print(p)


####################################################################################################
#Plot 5, 50 500 trees

#Plot the different trees of the Random Forest
brf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header=T)
prf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_rf_100_to_500tree_cv_full_chr.txt", header=T)

#Start by taking genes with best CV R2 > 0.3
#Then use the genes to subset the dataframe of random forest trees
brf <- subset(brf, CV_R2 > 0.5)
brf <- brf[,c(1,2)]
brf$Gene_ID <- as.character(brf$Gene_ID)

prf <- prf[,c(1,4:13)]
prf$gene_id <- as.character(prf$gene_id)

library(tidyverse)

rftrees <- inner_join(brf, prf, by = c("Gene_ID" = "gene_id"))
#I had to write out the dataframe to file so as to be able to use the "unique' function to take the unique genes in the dataframe
write.table(rftrees, file = "Z:/data/mesa_models/python_ml_models/rf_cv_gene_param.txt", quote = FALSE, sep = "\t", row.names = FALSE)
rftrees <- read.table(file = "Z:/data/mesa_models/python_ml_models/rf_cv_gene_param.txt", header=T)
rftrees$Gene_ID <- as.character(rftrees$Gene_ID)

rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/AFA_rf_5tree_cv_full_chr.txt", header=T, sep="\t")
rf_afa <- rf_afa[,c(1,3)]
rf_afa$Gene_ID <- as.character(rf_afa$Gene_ID)
rf5t <- rf_afa
rf5t$Gene_ID <- as.character(rf5t$Gene_ID)

rftrees <- inner_join(rf5t, rftrees, by = c("Gene_ID" = "Gene_ID"))

#I had to write out the dataframe to file so as to be able to use the "unique' function to take the unique genes in the dataframe
write.table(rftrees, file = "Z:/data/mesa_models/python_ml_models/rf_cv_gene_param.txt", quote = FALSE, sep = "\t", row.names = FALSE)
rftrees <- read.table(file = "Z:/data/mesa_models/python_ml_models/rf_cv_gene_param.txt", header=T)

rftrees <- rftrees[unique(rftrees$Gene_ID),]
names(rftrees) <- c('gene_id', '5', 'gene_name', '50','100','150','200','250','300','350','400','450','500')
rftrees <- rftrees[,c(2,4:13)]
rownames(rftrees) <- c(1:length(rftrees$`5`))

#Take out only 3 trees 5, 50, and 500
rf3 <- rftrees[,c(1,2,11)]

trees <- c(5,50,500)
trees <- as.factor(trees)

cl <- rainbow(length(rftrees$`5`))
plot(trees, rf3[1,], type = 'o', xlab = "Trees", ylab = "CV R2", main = "RF Tree", ylim = c(0.3,1), col=cl[1])

for (i in 2:length(rftrees$`50`)){
  lines(trees, rf3[i,], type = "o", col=cl[i+3])
}
barplot(as.matrix(rf3), legend.text = T, col = c('red','blue','black'), ylim = c(0.5,1))











############################################################################################
# Merge HIS best grids

"%&%" <- function(a,b) paste(a,b, sep = "")


#Merging all the chunks that were complete
algs <- c("rf", "svr", "knn")
pop <- "HIS"

#use this to merge all chunks into their singular chromosomes
#first check what each chunk file looks like
chk1rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/his_results/grid_split/HIS_best_grid_split_rf_cv_chr1_chunk1.txt", header=T)
chk2rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/his_results/grid_split/HIS_best_grid_split_rf_cv_chr1_chunk2.txt", header=T)
chk3rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/his_results/grid_split/HIS_best_grid_split_rf_cv_chr1_chunk3.txt", header=T)
chk4rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/his_results/grid_split/HIS_best_grid_split_rf_cv_chr1_chunk4.txt", header=T)
chk5rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/his_results/grid_split/HIS_best_grid_split_rf_cv_chr1_chunk5.txt", header=T)

allchk1 <- rbind(chk1rf, chk2rf, chk3rf, chk4rf, chk5rf)

for (alg in algs){
  for (no in 1:22){
    allchk <- NULL
    for (k in 1:5){
      allchk <- rbind(allchk, read.table(file = "Z:/data/mesa_models/python_ml_models/his_results/grid_split/" %&% pop %&% "_best_grid_split_"%&% alg %&% "_cv_chr" %&% no %&% "_chunk" %&% k %&% ".txt", header = T, sep = "\t", stringsAsFactors = F))
    }
    write.table(allchk, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/" %&% pop %&% "_best_grid_"%&% alg %&% "_chr" %&% no %&% "_full.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

bgrid_his_chr1_rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_chr1_full.txt", header = T)
bgrid_his_chr1_rf <- bgrid_his_chr1_rf[unique(bgrid_his_chr1_rf$Gene_Name),]

#This is the correct way I Should do merge all the chunks and the chromosomes into one single large file
for (alg in algs){
  allchr <- NULL
  for (no in 1:22){
    #allchr <- NULL
    for (k in 1:5){
      allchr <- rbind(allchr, read.table(file = "Z:/data/mesa_models/python_ml_models/his_results/grid_split/" %&% pop %&% "_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk" %&% k %&% ".txt", header = T, sep = "\t", stringsAsFactors = F))
    }
  }
  write.table(allchr, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/" %&% pop %&% "_best_grid_"%&% alg %&% "_all_chrom.txt", quote = FALSE, sep = "\t", row.names = FALSE)
}

bgrid_his_rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_rf_all_chrom.txt", header = T)
bgrid_his_knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_knn_all_chrom.txt", header = T)
bgrid_his_svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_svr_all_chrom.txt", header = T)

bgrid_his_knn <- bgrid_his_knn[unique(bgrid_his_knn$Gene_Name),]

#Because the total genes in all chroms was 9501 (should be 9623), I used the loop below to merge all chroms again
for (alg in algs){
  allchr <- NULL
  for (no in 1:22){
    allchr <- rbind(allchr, read.table(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_"%&% alg %&% "_chr" %&% no %&% "_full.txt", header = T, sep = "\t", stringsAsFactors = F))
  }
}

################################################################################################################################
# Merge HIS 2 METS

"%&%" <- function(a,b) paste(a,b, sep = "")

knn <- NULL
rf <- NULL
svr <- NULL
pop <- "HIS"

for (chrom in 1:22) {
  no <- as.character(chrom)
  knn <- rbind(knn, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_knn_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  rf <- rbind(rf, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_rf_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  svr <- rbind(svr, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_svr_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(knn, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_knn_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rf, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_rf_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_svr_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#############################################################################################################################
#Compare model Optimized performance on HIS 2 METS as against Elastic Net
#HIS 2METS

his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
his_2_mets$gene <- as.character(his_2_mets$gene) # 6512 genes
his_2_mets <- subset(his_2_mets, spearman > 0.1) #2146

knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_knn_cor_test_full_chr.txt", header = T)
knn <- knn[,c(1,9)]
knn$gene_id <- as.character(knn$gene_id) #8807 genes
knn <- subset(knn, knn$spearman_yobs_vs_ypred..d. > 0.1) #2239

rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_rf_cor_test_full_chr.txt", header = T)
rf <- rf[,c(1,9)]
rf$gene_id <- as.character(rf$gene_id) #8807
rf <- subset(rf, rf$spearman_yobs_vs_ypred..d. > 0.1) #2573

svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_svr_cor_test_full_chr.txt", header = T)
svr <- svr[,c(1,9)]
svr$gene_id <- as.character(svr$gene_id) #8807
svr <- subset(svr, svr$spearman_yobs_vs_ypred..d. > 0.1) #2416

#Take the overlapping genes between elastic net and each model
#RF
library(dplyr)
library(ggplot2)
library("ggpubr")
elrf <- inner_join(his_2_mets, rf, by = c("gene" = "gene_id")) #1168
elrf <- elrf[,c(2,3)]
names(elrf) <- c("elnet", "rf")

#Before plotting, take the difference of rf and elnet (rf - elnet) plot on y axis against elnet on x axis
rfdif <- elrf$rf - elrf$elnet

ggplot(elrf, aes(x=elnet, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("RF Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elrf, x = "elnet", y = "rf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "RF Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs (rf - elnet)
elrf <- cbind(elrf, rfdif)

ggplot(elrf, aes(x=elnet, y=rfdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("RF - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where rf minus elnet > 0. That is, how mnay genes where rf outperformed elnet. total genes = 6110
sum(elrf$rfdif > 0) # answer = 3064.Therefore elnet performed better only on 3046 genes
#after r >0.1 rf has 605 genes better

#Compare with non optimized and optimized rf
orf <- read.table(file="Z:/data/mesa_models/python_ml_models/results/HIS_2_METS_rf_cor_test_all_chr.txt", header = T, sep = "\t")
orf <- orf[,c(1,9)]
orf$gene_id <- as.character(orf$gene_id)

rforf <- inner_join(rf, orf, by = c("gene_id" = "gene_id"))
names(rforf) <- c("gene", "rf", "orf")

ggplot(rforf, aes(x=rf, y=orf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("Non Optimized RF") + xlab("Optimized RF") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(rforf, x = "rf", y = "orf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized RF", ylab = "Non Optimized RF", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

orfdif <- rforf$orf - rforf$rf
rforf <- cbind(rforf, orfdif) #I think raandom forest with fixed 100 trees is the most optimal

ggplot(rforf, aes(x=rf, y=orfdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("Non Optimized - Optimized RF") + xlab("Non Optimized") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#SVR
library(dplyr)
library(ggplot2)
library("ggpubr")
elsvr <- inner_join(his_2_mets, svr, by = c("gene" = "gene_id")) #1057
elsvr <- elsvr[,c(2,3)]
names(elsvr) <- c("elnet", "svr")
ggplot(elsvr, aes(x=elnet, y=svr)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elsvr, x = "elnet", y = "svr", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs SVR - Elnet
svrdif <- elsvr$svr - elsvr$elnet
elsvr <- cbind(elsvr, svrdif)

ggplot(elsvr, aes(x=elnet, y=svrdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where svr minus elnet > 0. That is, how many genes where svr outperformed elnet. total genes = 6110
sum(elsvr$svrdif > 0) # answer = 2908.Therefore elnet performed better only on 3202 genes
#after r > 0.1 487

#Compare with non optimized and optimized svr
osvr <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_svr_rbf_cor_test_all_chr.txt", header = T, sep = "\t")
osvr <- osvr[,c(1,9)]
osvr$gene_id <- as.character(osvr$gene_id)

#KNN
library(dplyr)
library(ggplot2)
library("ggpubr")
elknn <- inner_join(his_2_mets, knn, by = c("gene" = "gene_id")) #877
elknn <- elknn[,c(2,3)]
names(elknn) <- c("elnet", "knn")
ggplot(elknn, aes(x=elnet, y=knn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("KNN Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elknn, x = "elnet", y = "knn", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs knn - elnet
knndif <- elknn$knn - elknn$elnet
elknn <- cbind(elknn, knndif)
ggplot(elknn, aes(x=elnet, y=knndif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("KNN - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where knn minus elnet > 0. That is, how many genes where svr outperformed elnet. total genes = 6110
sum(elknn$knndif > 0) # answer = 2781.Therefore elnet performed better only on 3329 genes
# after r > 0.1 375-

#Compare with non optimized and optimized KNN
oknn <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_knn_cor_test_all_chr.txt", header = T, sep = "\t")
oknn <- oknn[,c(1,9)]
oknn$gene_id <- as.character(oknn$gene_id)

knnoknn <- inner_join(knn, oknn, by = c("gene_id" = "gene_id"))
names(knnoknn) <- c("gene", "knn", "oknn")

ggplot(knnoknn, aes(x=knn, y=oknn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("Non Optimized KNN") + xlab("Optimized KNN") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(knnoknn, x = "knn", y = "oknn", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized KNN", ylab = "Non Optimized KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")



###################################################################################################################################

his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
his_2_mets$gene <- as.character(his_2_mets$gene)

rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_rf_cor_test_full_chr.txt", header = T)
rf <- rf[,c(1,9)]
rf$gene_id <- as.character(rf$gene_id)



cau6_2 <- read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/CAU_best_grid_split_rf_cv_chr6_chunk2.txt", header = T)









##############################################################################################################################
# Merge CAU best grids

"%&%" <- function(a,b) paste(a,b, sep = "")


#Merging all the chunks that were complete
algs <- c("rf", "svr", "knn")
pop <- "CAU"

#use this to merge all chunks into their singular chromosomes
#first check what each chunk file looks like
# by the date 1 october 2019 9:57AM it is remaining 12 genes in chr 6 chunk 2 ie 82/104

chk1rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/CAU_best_grid_split_rf_cv_chr1_chunk1.txt", header=T)
chk2rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/CAU_best_grid_split_rf_cv_chr1_chunk2.txt", header=T)
chk3rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/CAU_best_grid_split_rf_cv_chr1_chunk3.txt", header=T)
chk4rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/CAU_best_grid_split_rf_cv_chr1_chunk4.txt", header=T)
chk5rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/CAU_best_grid_split_rf_cv_chr1_chunk5.txt", header=T)

allchk1 <- rbind(chk1rf, chk2rf, chk3rf, chk4rf, chk5rf)

for (alg in algs){
  for (no in 1:22){
    allchk <- NULL
    for (k in 1:5){
      allchk <- rbind(allchk, read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/" %&% pop %&% "_best_grid_split_"%&% alg %&% "_cv_chr" %&% no %&% "_chunk" %&% k %&% ".txt", header = T, sep = "\t", stringsAsFactors = F))
    }
    write.table(allchk, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/" %&% pop %&% "_best_grid_"%&% alg %&% "_chr" %&% no %&% "_full.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

bgrid_cau_chr1_rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_rf_chr1_full.txt", header = T)
bgrid_cau_chr1_rf <- bgrid_his_chr1_rf[unique(bgrid_his_chr1_rf$Gene_Name),]

#This is the correct way I Should do merge all the chunks and the chromosomes into one single large file
for (alg in algs){
  allchr <- NULL
  for (no in 1:22){
    #allchr <- NULL
    for (k in 1:5){
      allchr <- rbind(allchr, read.table(file = "Z:/data/mesa_models/python_ml_models/cau_results/grid_split/" %&% pop %&% "_best_grid_split_"%&% alg %&% "_cv_chr"%&% no %&% "_chunk" %&% k %&% ".txt", header = T, sep = "\t", stringsAsFactors = F))
    }
  }
  write.table(allchr, file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/" %&% pop %&% "_best_grid_"%&% alg %&% "_all_chrom.txt", quote = FALSE, sep = "\t", row.names = FALSE)
}

bgrid_cau_rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_rf_all_chrom.txt", header = T)
bgrid_cau_knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_knn_all_chrom.txt", header = T)
bgrid_cau_svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/CAU_best_grid_svr_all_chrom.txt", header = T)

bgrid_his_knn <- bgrid_his_knn[unique(bgrid_his_knn$Gene_Name),]

#Because the total genes in all chroms was 9479 (should be 9623), I used the loop below to merge all chroms again
for (alg in algs){
  allchr <- NULL
  for (no in 1:22){
    allchr <- rbind(allchr, read.table(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/HIS_best_grid_"%&% alg %&% "_chr" %&% no %&% "_full.txt", header = T, sep = "\t", stringsAsFactors = F))
  }
}

################################################################################################################################
# Merge CAU 2 METS

"%&%" <- function(a,b) paste(a,b, sep = "")

knn <- NULL
rf <- NULL
svr <- NULL
pop <- "CAU"

for (chrom in 1:22) {
  no <- as.character(chrom)
  knn <- rbind(knn, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_knn_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  rf <- rbind(rf, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_rf_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  svr <- rbind(svr, read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_svr_cor_test_chr" %&% chrom %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(knn, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_knn_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rf, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_rf_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr, file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_" %&% pop %&% "_2_METS_svr_cor_test_full_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


#############################################################################################################################
#Compare model Optimized performance on HIS 2 METS as against Elastic Net
#HIS 2METS

his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
his_2_mets$gene <- as.character(his_2_mets$gene)

knn <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_knn_cor_test_full_chr.txt", header = T)
knn <- knn[,c(1,9)]
knn$gene_id <- as.character(knn$gene_id)

rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_rf_cor_test_full_chr.txt", header = T)
rf <- rf[,c(1,9)]
rf$gene_id <- as.character(rf$gene_id)

svr <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_svr_cor_test_full_chr.txt", header = T)
svr <- svr[,c(1,9)]
svr$gene_id <- as.character(svr$gene_id)

#Take the overlapping genes between elastic net and each model
#RF
library(dplyr)
library(ggplot2)
library("ggpubr")
elrf <- inner_join(his_2_mets, rf, by = c("gene" = "gene_id"))
elrf <- elrf[,c(2,3)]
names(elrf) <- c("elnet", "rf")

#Before plotting, take the difference of rf and elnet (rf - elnet) plot on y axis against elnet on x axis
rfdif <- elrf$rf - elrf$elnet

ggplot(elrf, aes(x=elnet, y=rf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("RF Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elrf, x = "elnet", y = "rf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "RF Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs (rf - elnet)
elrf <- cbind(elrf, rfdif)

ggplot(elrf, aes(x=elnet, y=rfdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("RF - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where rf minus elnet > 0. That is, how mnay genes where rf outperformed elnet. total genes = 6110
sum(elrf$rfdif > 0) # answer = 3064.Therefore elnet performed better only on 3046 genes

#Compare with non optimized and optimized rf
orf <- read.table(file="Z:/data/mesa_models/python_ml_models/results/HIS_2_METS_rf_cor_test_all_chr.txt", header = T, sep = "\t")
orf <- orf[,c(1,9)]
orf$gene_id <- as.character(orf$gene_id)

rforf <- inner_join(rf, orf, by = c("gene_id" = "gene_id"))
names(rforf) <- c("gene", "rf", "orf")

ggplot(rforf, aes(x=rf, y=orf)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("Non Optimized RF") + xlab("Optimized RF") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(rforf, x = "rf", y = "orf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized RF", ylab = "Non Optimized RF", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

orfdif <- rforf$orf - rforf$rf
rforf <- cbind(rforf, orfdif) #I think raandom forest with fixed 100 trees is the most optimal

ggplot(rforf, aes(x=rf, y=orfdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("Non Optimized - Optimized RF") + xlab("Non Optimized") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#SVR
library(dplyr)
library(ggplot2)
library("ggpubr")
elsvr <- inner_join(his_2_mets, svr, by = c("gene" = "gene_id"))
elsvr <- elsvr[,c(2,3)]
names(elsvr) <- c("elnet", "svr")
ggplot(elsvr, aes(x=elnet, y=svr)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elsvr, x = "elnet", y = "svr", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs SVR - Elnet
svrdif <- elsvr$svr - elsvr$elnet
elsvr <- cbind(elsvr, svrdif)

ggplot(elsvr, aes(x=elnet, y=svrdif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where svr minus elnet > 0. That is, how many genes where svr outperformed elnet. total genes = 6110
sum(elsvr$svrdif > 0) # answer = 2908.Therefore elnet performed better only on 3202 genes

#Compare with non optimized and optimized svr
osvr <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_svr_rbf_cor_test_all_chr.txt", header = T, sep = "\t")
osvr <- osvr[,c(1,9)]
osvr$gene_id <- as.character(osvr$gene_id)

#KNN
library(dplyr)
library(ggplot2)
library("ggpubr")
elknn <- inner_join(his_2_mets, knn, by = c("gene" = "gene_id"))
elknn <- elknn[,c(2,3)]
names(elknn) <- c("elnet", "knn")
ggplot(elknn, aes(x=elnet, y=knn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("KNN Optimized") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(elknn, x = "elnet", y = "knn", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#Plot elnet vs knn - elnet
knndif <- elknn$knn - elknn$elnet
elknn <- cbind(elknn, knndif)
ggplot(elknn, aes(x=elnet, y=knndif)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("KNN - Elastic Net") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=0,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

#Count how many genes where knn minus elnet > 0. That is, how many genes where svr outperformed elnet. total genes = 6110
sum(elknn$knndif > 0) # answer = 2781.Therefore elnet performed better only on 3329 genes

#Compare with non optimized and optimized KNN
oknn <- read.table(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_knn_cor_test_all_chr.txt", header = T, sep = "\t")
oknn <- oknn[,c(1,9)]
oknn$gene_id <- as.character(oknn$gene_id)

knnoknn <- inner_join(knn, oknn, by = c("gene_id" = "gene_id"))
names(knnoknn) <- c("gene", "knn", "oknn")

ggplot(knnoknn, aes(x=knn, y=oknn)) + ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("Non Optimized KNN") + xlab("Optimized KNN") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-1,1)) + ylim(c(-1,1))

ggscatter(knnoknn, x = "knn", y = "oknn", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Optimized KNN", ylab = "Non Optimized KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")



###################################################################################################################################

his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
his_2_mets$gene <- as.character(his_2_mets$gene)

rf <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_HIS_2_METS_rf_cor_test_full_chr.txt", header = T)
rf <- rf[,c(1,9)]
rf$gene_id <- as.character(rf$gene_id)














###############################################################################################################################
#filter the ML models with CV R2 > 0.01
library(dplyr)
library(ggplot2)
library("ggpubr")
#AFA
#elnet_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE) #elnet
elnet_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
elnet_afa_2_mets$gene <- as.character(elnet_afa_2_mets$gene)

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

elnet_afa_2_mets_0.1 <- subset(elnet_afa_2_mets, spearman > 0.1)
filt_rf_afa_2_mets_0.1 <- subset(filt_rf_afa_2_mets, rf_spearman > 0.1)
filt_elnet_rf_afa_2_mets_0.1 <- inner_join(elnet_afa_2_mets_0.1, filt_rf_afa_2_mets_0.1, by = c("gene" = "gene"))
#plot, NOTE all models here are optimized
ggplot(filt_elnet_rf_afa_2_mets_0.1, aes(x=spearman, y=rf_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("Random Forest") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1)) + 
  theme_bw(20)

ggscatter(filt_elnet_rf_afa_2_mets_0.1, x = "spearman", y = "rf_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "Random Forest",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")+
          theme_classic2(20)
  

filt_elnet_rf_afa_2_mets <- inner_join(elnet_afa_2_mets, filt_rf_afa_2_mets, by = c("gene" = "gene"))
sub_elnet_rf_afa_2_mets <- subset(filt_elnet_rf_afa_2_mets, spearman > 0.1 | rf_spearman > 0.1) #removes where both are < 0.1
#suband_elnet_rf_afa_2_mets <- subset(filt_elnet_rf_afa_2_mets, spearman > 0.1 & rf_spearman > 0.1) # where both are > 0.1
#Plot
ggplot(sub_elnet_rf_afa_2_mets, aes(x=spearman, y=rf_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("Random Forest") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.4,1)) + ylim(c(-0.4,1)) + theme_bw()

#Then find the genes where EN or RF is positive and the other negative 
#Do FUMA on top 10 of theses genes
rfpos <- subset(sub_elnet_rf_afa_2_mets, spearman < 0)
enpos <- subset(sub_elnet_rf_afa_2_mets, rf_spearman < 0)

##########
#get the genes that are not in the other. rf genes not in en. also en genes not in rf
#Do FUMA on top 10 of these genes
rfonly <- left_join(filt_rf_afa_2_mets_0.1, elnet_afa_2_mets_0.1, by = c("gene" = "gene"))
rfonly <- rfonly[is.na(rfonly$spearman),]
rfonly <- anti_join(filt_rf_afa_2_mets_0.1, elnet_afa_2_mets_0.1, by = c("gene" = "gene")) #Better!
rfonly_den <- data.frame(spearman=rfonly$rf_spearman, prediction=rep("RF Only", length(rfonly$rf_spearman)))


enonly <- left_join(elnet_afa_2_mets_0.1, filt_rf_afa_2_mets_0.1,by = c("gene" = "gene"))
enonly <- enonly[is.na(enonly$rf_spearman),]
enonly <- anti_join(elnet_afa_2_mets_0.1, filt_rf_afa_2_mets_0.1,by = c("gene" = "gene"))
enonly_den <- data.frame(spearman=enonly$spearman, prediction=rep("EN Only", length(enonly$spearman)))


enplusrf_rf <- data.frame(spearman=filt_elnet_rf_afa_2_mets_0.1$rf_spearman, 
                          prediction=rep("EN+RF RF Only",length(filt_elnet_rf_afa_2_mets_0.1$rf_spearman)))

enplusrf_en <- data.frame(spearman=filt_elnet_rf_afa_2_mets_0.1$spearman, 
                          prediction=rep("EN+RF EN Only",length(filt_elnet_rf_afa_2_mets_0.1$spearman)))


den_afa_2_mets <- rbind(rfonly_den, enonly_den, enplusrf_rf, enplusrf_en)

enonly_e <- anti_join(enonly,filt_elnet_rf_afa_2_mets_0.1,by = c("gene" = "gene"))
#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("Random Forest", "Elastic NET", "Intersect"), no=c(239, 668, 928)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (AFA 2 METS)") + geom_text(aes(label=no), vjust=3, color="white", size=10) +
  scale_fill_manual(values = c("blue","violet","red")) + 
  theme_minimal(20) # scale_fill_manual is how to give my own color specification

#Density plots R > 0.1
ggplot(enonly, aes(x = spearman)) + scale_x_continuous(name = "Spearman Correlation > 0.1") +
  scale_y_continuous(name = "Density") + ggtitle("Density plot of Genes only in Elastic Net (AFA 2 METS)") +
  geom_density(fill = "red", colour = "red", alpha=0.6) + theme_bw() +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")#EN only genes

ggplot(rfonly, aes(x = rf_spearman)) + scale_x_continuous(name = "Spearman Correlation > 0.1") +
  scale_y_continuous(name = "Density") + ggtitle("Density plot of Genes only in Random Forest (AFA 2 METS)") +
  geom_density(fill = "blue", colour = "blue", alpha=0.6) + theme_bw() +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")#RF only genes

ggplot(filt_elnet_rf_afa_2_mets_0.1, aes(x = rf_spearman)) + scale_x_continuous(name = "Spearman Correlation > 0.1") +
  scale_y_continuous(name = "Density") + ggtitle("Density plot of Genes in Intersect of RF and EN for Random Forest (AFA 2 METS)") +
  geom_density(fill = "green", colour = "green", alpha=0.6) + theme_bw() +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed") #intersect of RF and EN

ggplot(filt_elnet_rf_afa_2_mets_0.1, aes(x = spearman)) + scale_x_continuous(name = "Spearman Correlation > 0.1") +
  scale_y_continuous(name = "Density") + ggtitle("Density plot of Genes in Intersect of RF and EN for Elastic Net (AFA 2 METS)") +
  geom_density(fill = "green4", colour = "green4", alpha=0.6) + theme_bw() +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed") #intersect of RF and EN

#find the mean of each ML group
library(plyr)
mu <- ddply(den_afa_2_mets, "prediction", summarise, grp.mean=mean(spearman))

#Combine all density plot into one
ggplot(den_afa_2_mets, aes(x = spearman, color=prediction)) + scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Genes in Intersect of RF and EN for Elastic Net (AFA 2 METS)") +
  geom_density()+ theme_bw() +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")

# lwd = line thickness
ggplot(den_afa_2_mets, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") +
  geom_density()+ theme_classic(20) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed") + 
  scale_color_manual(values = c("red","blue","orange","violet"))


# make density plot for EN intersect RF AFA 2 METS
int_enrf_en <- data.frame(spearman=filt_elnet_rf_afa_2_mets_0.1$spearman, 
                          prediction=rep("EN",length(filt_elnet_rf_afa_2_mets_0.1$spearman)))

int_enrf_rf <- data.frame(spearman=filt_elnet_rf_afa_2_mets_0.1$rf_spearman, 
                          prediction=rep("RF",length(filt_elnet_rf_afa_2_mets_0.1$rf_spearman)))

den_int_enrf <- rbind(int_enrf_en, int_enrf_rf)
#find the mean of each ML group
library(plyr)
mu <- ddply(den_int_enrf, "prediction", summarise, grp.mean=mean(spearman))

# lwd = line thickness
ggplot(den_int_enrf, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Intersect Genes for RF and EN (AFA 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")

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
names(filt_svr_afa_2_mets) <- c("gene", "svr_spearman")

elnet_afa_2_mets_0.1 <- subset(elnet_afa_2_mets, spearman > 0.1)
filt_svr_afa_2_mets_0.1 <- subset(filt_svr_afa_2_mets, svr_spearman > 0.1)
filt_elnet_svr_afa_2_mets_0.1 <- inner_join(elnet_afa_2_mets_0.1, filt_svr_afa_2_mets_0.1, by = c("gene" = "gene"))
#plot, NOTE all models here are optimized
ggplot(filt_elnet_svr_afa_2_mets_0.1, aes(x=spearman, y=svr_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("SVR") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_svr_afa_2_mets_0.1, x = "spearman", y = "svr_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")+
          theme_classic2(20)

#intersect before filtering off R > 0.1 for each of ML and EN
filt_elnet_svr_afa_2_mets <- inner_join(elnet_afa_2_mets, filt_svr_afa_2_mets, by = c("gene" = "gene"))
sub_elnet_svr_afa_2_mets <- subset(filt_elnet_svr_afa_2_mets, spearman > 0.1 | svr_spearman > 0.1) #then filter off where both are < 0.1

#Plot
ggplot(sub_elnet_svr_afa_2_mets, aes(x=spearman, y=svr_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("SVR") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.4,1)) + ylim(c(-0.4,1)) + theme_bw()
#Then find the genes where EN or SVR is positive and the other negative 
#Do FUMA on top 10 of theses genes
svrpos <- subset(sub_elnet_svr_afa_2_mets, spearman < 0)
enpos_s <- subset(sub_elnet_svr_afa_2_mets, svr_spearman < 0)

##########
#get the genes that are not in the other. svr genes not in en. also en genes not in svr
#Do FUMA on top 10 of these genes
svronly <- anti_join(filt_svr_afa_2_mets_0.1, elnet_afa_2_mets_0.1, by = c("gene" = "gene")) #Better!

enonly_s <- anti_join(elnet_afa_2_mets_0.1, filt_svr_afa_2_mets_0.1,by = c("gene" = "gene"))

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("SVR", "Elastic NET", "Intersect"), no=c(270, 904, 691)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (AFA 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here


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
names(filt_knn_afa_2_mets) <- c("gene", "knn_spearman")

elnet_afa_2_mets_0.1 <- subset(elnet_afa_2_mets, spearman > 0.1)
filt_knn_afa_2_mets_0.1 <- subset(filt_knn_afa_2_mets, knn_spearman > 0.1)
filt_elnet_knn_afa_2_mets_0.1 <- inner_join(elnet_afa_2_mets_0.1, filt_knn_afa_2_mets_0.1, by = c("gene" = "gene"))
#plot Intersect where both ML R > 0.1, NOTE all models here are optimized
ggplot(filt_elnet_knn_afa_2_mets_0.1, aes(x=spearman, y=knn_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_knn_afa_2_mets_0.1, x = "spearman", y = "knn_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")+
          theme_classic2(20)

#intersect before filtering off R > 0.1 for each of ML and EN
filt_elnet_knn_afa_2_mets <- inner_join(elnet_afa_2_mets, filt_knn_afa_2_mets, by = c("gene" = "gene"))
sub_elnet_knn_afa_2_mets <- subset(filt_elnet_knn_afa_2_mets, spearman > 0.1 | knn_spearman > 0.1) #then filter off where both are < 0.1

#Plot
ggplot(sub_elnet_knn_afa_2_mets, aes(x=spearman, y=knn_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.4,1)) + ylim(c(-0.4,1)) + theme_bw()
#Then find the genes where EN or KNN is positive and the other negative 
#Do FUMA on top 10 of theses genes
knnpos <- subset(sub_elnet_knn_afa_2_mets, spearman < 0)
enpos_k <- subset(sub_elnet_knn_afa_2_mets, knn_spearman < 0)

##########
#get the genes that are not in the other. svr genes not in en. also en genes not in svr
#Do FUMA on top 10 of these genes
knnonly <- anti_join(filt_knn_afa_2_mets_0.1, elnet_afa_2_mets_0.1, by = c("gene" = "gene")) #Better!

enonly_k <- anti_join(elnet_afa_2_mets_0.1, filt_knn_afa_2_mets_0.1,by = c("gene" = "gene"))

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("KNN", "Elastic NET", "Intersect"), no=c(264, 1033, 560)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (AFA 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here

#Make bar plot with count of genes in only ML, EN+ML, and EN for genes R > 0.1
gcount <- data.frame(ml=c("EN-RF", "EN+RF", "RF-EN", "EN-SVR", "EN+SVR", "SVR-EN",
                          "EN-KNN", "EN+KNN", "KNN-EN"), no=c(668,928,239,904,691,270,1033,560,264)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (AFA 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here

#Make violin plot for all 3 population mesa to mets with EN for R > 0.1
elnet_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
elnet_afa_2_mets$gene <- as.character(elnet_afa_2_mets$gene)
en_afa_2mets <- data.frame(spearman=elnet_afa_2_mets_0.1$spearman, prediction=rep("AFA_2_METS", length(elnet_afa_2_mets_0.1$spearman)))

elnet_his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
elnet_his_2_mets$gene <- as.character(elnet_his_2_mets$gene)
elnet_his_2_mets_0.1 <- subset(elnet_his_2_mets, spearman > 0.1)
en_his_2mets <- data.frame(spearman=elnet_his_2_mets_0.1$spearman, prediction=rep("HIS_2_METS", length(elnet_his_2_mets_0.1$spearman)))

elnet_cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
elnet_cau_2_mets$gene <- as.character(elnet_cau_2_mets$gene)
elnet_cau_2_mets_0.1 <- subset(elnet_cau_2_mets, spearman > 0.1)
en_cau_2mets <- data.frame(spearman=elnet_cau_2_mets_0.1$spearman, prediction=rep("CAU_2_METS", length(elnet_cau_2_mets_0.1$spearman)))

en_mesa_2_mets <- rbind(en_afa_2mets, en_his_2mets, en_cau_2mets)

ggplot(en_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, 
                           fill=prediction)) + 
  geom_violin(trim = T) + geom_boxplot(width=0.2, color="black") +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.2), limits=c(0, 1.0))


#Make violin plot for all 3 population mesa to mets with EN without doing R > 0.1
elnet_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
elnet_afa_2_mets$gene <- as.character(elnet_afa_2_mets$gene)
en_afa_2mets <- data.frame(spearman=elnet_afa_2_mets$spearman, 
                           prediction=rep("AFA_2_METS", length(elnet_afa_2_mets$spearman)))

elnet_his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
elnet_his_2_mets$gene <- as.character(elnet_his_2_mets$gene)
en_his_2mets <- data.frame(spearman=elnet_his_2_mets$spearman, 
                           prediction=rep("HIS_2_METS", length(elnet_his_2_mets$spearman)))

elnet_cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
elnet_cau_2_mets$gene <- as.character(elnet_cau_2_mets$gene)
en_cau_2mets <- data.frame(spearman=elnet_cau_2_mets$spearman, 
                           prediction=rep("CAU_2_METS", length(elnet_cau_2_mets$spearman)))

en_mesa_2_mets <- rbind(en_afa_2mets, en_his_2mets, en_cau_2mets)
# Function to produce summary statistics (mean and +/- sd)

ggplot(en_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction)) + 
  geom_violin(trim = F) + geom_boxplot(width=0.2, color="black") + theme_minimal() + 
  ggtitle("Elastic Net (MESA to METS)") +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw() +
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.2), limits=c(0, 1.0))



#Make density plot for all 4 ML models with ML R > 0.1 AFA 2 METS
rf_0.1 <- data.frame(spearman=filt_rf_afa_2_mets_0.1$rf_spearman, 
                     prediction=rep("RF", length(filt_rf_afa_2_mets_0.1$rf_spearman)))

svr_0.1 <- data.frame(spearman=filt_svr_afa_2_mets_0.1$svr_spearman, 
                     prediction=rep("SVR", length(filt_svr_afa_2_mets_0.1$svr_spearman)))

knn_0.1 <- data.frame(spearman=filt_knn_afa_2_mets_0.1$knn_spearman, 
                     prediction=rep("KNN", length(filt_knn_afa_2_mets_0.1$knn_spearman)))

en_0.1 <- data.frame(spearman=elnet_afa_2_mets_0.1$spearman, 
                     prediction=rep("EN", length(elnet_afa_2_mets_0.1$spearman)))
den_4_ml <- rbind(rf_0.1, svr_0.1, knn_0.1, en_0.1)


#find the mean of each ML group
library(plyr)
mu <- ddply(den_4_ml, "prediction", summarise, grp.mean=mean(spearman))

# density plot for all ML in afa 2 mets
# lwd = line thickness
ggplot(den_4_ml, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Genes for All Models (AFA 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")








#HIS

###############################################################################################################################
#filter the ML models with CV R2 > 0.01
library(dplyr)
library(ggplot2)
library("ggpubr")
#HIS
#elnet_afa <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE) #elnet
elnet_his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
elnet_his_2_mets$gene <- as.character(elnet_his_2_mets$gene)

#RF
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
names(filt_rf_his_2_mets) <- c("gene", "rf_spearman")

elnet_his_2_mets_0.1 <- subset(elnet_his_2_mets, spearman > 0.1)
filt_rf_his_2_mets_0.1 <- subset(filt_rf_his_2_mets, rf_spearman > 0.1)
filt_elnet_rf_his_2_mets_0.1 <- inner_join(elnet_his_2_mets_0.1, filt_rf_his_2_mets_0.1, by = c("gene" = "gene"))
#plot, NOTE all models here are optimized
ggplot(filt_elnet_rf_his_2_mets_0.1, aes(x=spearman, y=rf_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("Random Forest") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_rf_his_2_mets_0.1, x = "spearman", y = "rf_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "Random Forest", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

filt_elnet_rf_his_2_mets <- inner_join(elnet_his_2_mets, filt_rf_his_2_mets, by = c("gene" = "gene"))
sub_elnet_rf_his_2_mets <- subset(filt_elnet_rf_his_2_mets, spearman > 0.1 | rf_spearman > 0.1) #removes where both are < 0.1
#suband_elnet_rf_afa_2_mets <- subset(filt_elnet_rf_afa_2_mets, spearman > 0.1 & rf_spearman > 0.1) # where both are > 0.1
#Plot
ggplot(sub_elnet_rf_his_2_mets, aes(x=spearman, y=rf_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)") + 
  ylab("Random Forest") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.4,1)) + ylim(c(-0.4,1)) + theme_bw()

#Then find the genes where EN or RF is positive and the other negative 
#Do FUMA on top 10 of theses genes
rfpos_his <- subset(sub_elnet_rf_his_2_mets, spearman < 0)
enpos_his <- subset(sub_elnet_rf_his_2_mets, rf_spearman < 0)

##########
#get the genes that are not in the other. rf genes not in en. also en genes not in rf
#Do FUMA on top 10 of these genes
rfonly_his <- anti_join(filt_rf_his_2_mets_0.1, elnet_his_2_mets_0.1, by = c("gene" = "gene")) #Better!
rfonly_den_his <- data.frame(spearman=rfonly_his$rf_spearman, prediction=rep("RF Only", length(rfonly_his$rf_spearman)))

enonly_his <- anti_join(elnet_his_2_mets_0.1, filt_rf_his_2_mets_0.1,by = c("gene" = "gene"))
enonly_den_his <- data.frame(spearman=enonly_his$spearman, prediction=rep("EN Only", length(enonly_his$spearman)))


enplusrf_rf_his <- data.frame(spearman=filt_elnet_rf_his_2_mets_0.1$rf_spearman, 
                          prediction=rep("EN+RF RF Only",length(filt_elnet_rf_his_2_mets_0.1$rf_spearman)))

enplusrf_en_his <- data.frame(spearman=filt_elnet_rf_his_2_mets_0.1$spearman, 
                          prediction=rep("EN+RF EN Only",length(filt_elnet_rf_his_2_mets_0.1$spearman)))


den_his_2_mets <- rbind(rfonly_den_his, enonly_den_his, enplusrf_rf_his, enplusrf_en_his)

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("Random Forest", "Elastic NET", "Intersect"), no=c(355, 1356, 790)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (HIS 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine svr and knn here

#make density plot
#find the mean of each ML group
library(plyr)
mu <- ddply(den_his_2_mets, "prediction", summarise, grp.mean=mean(spearman))

#Combine all density plot into one
# lwd = line thickness
ggplot(den_his_2_mets, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Genes for RF and EN (HIS 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")


# make density plot for EN intersect RF HIS 2 METS
int_enrf_en_his <- data.frame(spearman=filt_elnet_rf_his_2_mets_0.1$spearman, 
                          prediction=rep("EN",length(filt_elnet_rf_his_2_mets_0.1$spearman)))

int_enrf_rf_his <- data.frame(spearman=filt_elnet_rf_his_2_mets_0.1$rf_spearman, 
                          prediction=rep("RF",length(filt_elnet_rf_his_2_mets_0.1$rf_spearman)))

den_int_enrf_his <- rbind(int_enrf_en_his, int_enrf_rf_his)
#find the mean of each ML group
library(plyr)
mu <- ddply(den_int_enrf_his, "prediction", summarise, grp.mean=mean(spearman))

# lwd = line thickness
ggplot(den_int_enrf_his, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Intersect Genes for RF and EN (HIS 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")

#SVR
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
names(filt_svr_his_2_mets) <- c("gene", "svr_spearman")

elnet_his_2_mets_0.1 <- subset(elnet_his_2_mets, spearman > 0.1)
filt_svr_his_2_mets_0.1 <- subset(filt_svr_his_2_mets, svr_spearman > 0.1)
filt_elnet_svr_his_2_mets_0.1 <- inner_join(elnet_his_2_mets_0.1, filt_svr_his_2_mets_0.1, by = c("gene" = "gene"))
#plot, NOTE all models here are optimized
ggplot(filt_elnet_svr_his_2_mets_0.1, aes(x=spearman, y=svr_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_svr_his_2_mets_0.1, x = "spearman", y = "svr_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#intersect before filtering off R > 0.1 for each of ML and EN
filt_elnet_svr_his_2_mets <- inner_join(elnet_his_2_mets, filt_svr_his_2_mets, by = c("gene" = "gene"))
sub_elnet_svr_his_2_mets <- subset(filt_elnet_svr_his_2_mets, spearman > 0.1 | svr_spearman > 0.1) #then filter off where both are < 0.1

#Plot
ggplot(sub_elnet_svr_his_2_mets, aes(x=spearman, y=svr_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("SVR") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.4,1)) + ylim(c(-0.4,1)) + theme_bw()
#Then find the genes where EN or SVR is positive and the other negative 
#Do FUMA on top 10 of theses genes
svrpos_his <- subset(sub_elnet_svr_his_2_mets, spearman < 0)
enpos_s_his <- subset(sub_elnet_svr_his_2_mets, svr_spearman < 0)

##########
#get the genes that are not in the other. svr genes not in en. also en genes not in svr
#Do FUMA on top 10 of these genes
svronly_his <- anti_join(filt_svr_his_2_mets_0.1, elnet_his_2_mets_0.1, by = c("gene" = "gene")) #Better!

enonly_s_his <- anti_join(elnet_his_2_mets_0.1, filt_svr_his_2_mets_0.1,by = c("gene" = "gene"))

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("SVR", "Elastic NET", "Intersect"), no=c(347, 1442, 704)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (HIS 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here


#KNN
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
names(filt_knn_his_2_mets) <- c("gene", "knn_spearman")

elnet_his_2_mets_0.1 <- subset(elnet_his_2_mets, spearman > 0.1)
filt_knn_his_2_mets_0.1 <- subset(filt_knn_his_2_mets, knn_spearman > 0.1)
filt_elnet_knn_his_2_mets_0.1 <- inner_join(elnet_his_2_mets_0.1, filt_knn_his_2_mets_0.1, by = c("gene" = "gene"))
#plot Intersect where both ML R > 0.1, NOTE all models here are optimized
ggplot(filt_elnet_knn_his_2_mets_0.1, aes(x=spearman, y=knn_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_knn_his_2_mets_0.1, x = "spearman", y = "knn_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#intersect before filtering off R > 0.1 for each of ML and EN
filt_elnet_knn_his_2_mets <- inner_join(elnet_his_2_mets, filt_knn_his_2_mets, by = c("gene" = "gene"))
sub_elnet_knn_his_2_mets <- subset(filt_elnet_knn_his_2_mets, spearman > 0.1 | knn_spearman > 0.1) #then filter off where both are < 0.1

#Plot
ggplot(sub_elnet_knn_his_2_mets, aes(x=spearman, y=knn_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (HIS to METS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.5,1)) + ylim(c(-0.5,1)) + theme_bw()
#Then find the genes where EN or KNN is positive and the other negative 
#Do FUMA on top 10 of theses genes
knnpos_his <- subset(sub_elnet_knn_his_2_mets, spearman < 0)
enpos_k_his <- subset(sub_elnet_knn_his_2_mets, knn_spearman < 0)

##########
#get the genes that are not in the other. svr genes not in en. also en genes not in svr
#Do FUMA on top 10 of these genes
knnonly_his <- anti_join(filt_knn_his_2_mets_0.1, elnet_his_2_mets_0.1, by = c("gene" = "gene")) #Better!

enonly_k_his <- anti_join(elnet_his_2_mets_0.1, filt_knn_his_2_mets_0.1,by = c("gene" = "gene"))

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("KNN", "Elastic NET", "Intersect"), no=c(268, 1628, 518)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (HIS 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here

#Make bar plot with count of genes in only ML, EN+ML, and EN for genes R > 0.1
gcount <- data.frame(ml=c("EN-RF", "EN+RF", "RF-EN", "EN-SVR", "EN+SVR", "SVR-EN",
                          "EN-KNN", "EN+KNN", "KNN-EN"), no=c(1356,790,355,1442,704,347,1628,518,268)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (HIS 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here

#Make violin plot for all 3 population mesa to mets with EN
elnet_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
elnet_afa_2_mets$gene <- as.character(elnet_afa_2_mets$gene)
en_afa_2mets <- data.frame(spearman=elnet_afa_2_mets_0.1$spearman, prediction=rep("AFA_2_METS", length(elnet_afa_2_mets_0.1$spearman)))

elnet_his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
elnet_his_2_mets$gene <- as.character(elnet_his_2_mets$gene)
elnet_his_2_mets_0.1 <- subset(elnet_his_2_mets, spearman > 0.1)
en_his_2mets <- data.frame(spearman=elnet_his_2_mets_0.1$spearman, prediction=rep("HIS_2_METS", length(elnet_his_2_mets_0.1$spearman)))

elnet_cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
elnet_cau_2_mets$gene <- as.character(elnet_cau_2_mets$gene)
elnet_cau_2_mets_0.1 <- subset(elnet_cau_2_mets, spearman > 0.1)
en_cau_2mets <- data.frame(spearman=elnet_cau_2_mets_0.1$spearman, prediction=rep("CAU_2_METS", length(elnet_cau_2_mets_0.1$spearman)))

en_mesa_2_mets <- rbind(en_afa_2mets, en_his_2mets, en_cau_2mets)
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

library(Hmisc)
ggplot(en_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, 
                           fill=prediction)) + 
  geom_violin(trim = F) + geom_boxplot(width=0.2, color="black") + theme_minimal() + 
  ggtitle("Elastic Net (MESA to METS)") +
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.2), limits=c(0, 1.0))#+
#stat_summary(fun.y=mean, geom = "point", size=3, color="white")

#Make density plot for all 4 ML models with ML R > 0.1 HIS 2 METS
rf_0.1_his <- data.frame(spearman=filt_rf_his_2_mets_0.1$rf_spearman, 
                     prediction=rep("RF", length(filt_rf_his_2_mets_0.1$rf_spearman)))

svr_0.1_his <- data.frame(spearman=filt_svr_his_2_mets_0.1$svr_spearman, 
                      prediction=rep("SVR", length(filt_svr_his_2_mets_0.1$svr_spearman)))

knn_0.1_his <- data.frame(spearman=filt_knn_his_2_mets_0.1$knn_spearman, 
                      prediction=rep("KNN", length(filt_knn_his_2_mets_0.1$knn_spearman)))

en_0.1_his <- data.frame(spearman=elnet_his_2_mets_0.1$spearman, 
                     prediction=rep("EN", length(elnet_his_2_mets_0.1$spearman)))
den_4_ml_his <- rbind(rf_0.1_his, svr_0.1_his, knn_0.1_his, en_0.1_his)


#find the mean of each ML group
library(plyr)
mu <- ddply(den_4_ml_his, "prediction", summarise, grp.mean=mean(spearman))

# density plot for all ML in afa 2 mets
# lwd = line thickness
ggplot(den_4_ml_his, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Genes for All Models (HIS 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")















#CAU

###############################################################################################################################
#filter the ML models with CV R2 > 0.01
library(dplyr)
library(ggplot2)
library("ggpubr")
#CAU
elnet_cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
elnet_cau_2_mets$gene <- as.character(elnet_cau_2_mets$gene)

#RF
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
names(filt_rf_cau_2_mets) <- c("gene", "rf_spearman")

elnet_cau_2_mets_0.1 <- subset(elnet_cau_2_mets, spearman > 0.1)
filt_rf_cau_2_mets_0.1 <- subset(filt_rf_cau_2_mets, rf_spearman > 0.1)
filt_elnet_rf_cau_2_mets_0.1 <- inner_join(elnet_cau_2_mets_0.1, filt_rf_cau_2_mets_0.1, by = c("gene" = "gene"))
#plot, NOTE all models here are optimized
ggplot(filt_elnet_rf_cau_2_mets_0.1, aes(x=spearman, y=rf_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)") + 
  ylab("Random Forest") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_rf_cau_2_mets_0.1, x = "spearman", y = "rf_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "Random Forest", title = "Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

filt_elnet_rf_cau_2_mets <- inner_join(elnet_cau_2_mets, filt_rf_cau_2_mets, by = c("gene" = "gene"))
sub_elnet_rf_cau_2_mets <- subset(filt_elnet_rf_cau_2_mets, spearman > 0.1 | rf_spearman > 0.1) #removes where both are < 0.1
#suband_elnet_rf_afa_2_mets <- subset(filt_elnet_rf_afa_2_mets, spearman > 0.1 & rf_spearman > 0.1) # where both are > 0.1
#Plot
ggplot(sub_elnet_rf_cau_2_mets, aes(x=spearman, y=rf_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)") + 
  ylab("Random Forest") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.4,1)) + ylim(c(-0.4,1)) + theme_bw()

#Then find the genes where EN or RF is positive and the other negative 
#Do FUMA on top 10 of theses genes
rfpos_cau <- subset(sub_elnet_rf_cau_2_mets, spearman < 0)
enpos_cau <- subset(sub_elnet_rf_cau_2_mets, rf_spearman < 0)

##########
#get the genes that are not in the other. rf genes not in en. also en genes not in rf
#Do FUMA on top 10 of these genes
rfonly_cau <- anti_join(filt_rf_cau_2_mets_0.1, elnet_cau_2_mets_0.1, by = c("gene" = "gene")) #Better!
rfonly_den_cau <- data.frame(spearman=rfonly_cau$rf_spearman, prediction=rep("RF Only", length(rfonly_cau$rf_spearman)))

enonly_cau <- anti_join(elnet_cau_2_mets_0.1, filt_rf_cau_2_mets_0.1,by = c("gene" = "gene"))
enonly_den_cau <- data.frame(spearman=enonly_cau$spearman, prediction=rep("EN Only", length(enonly_cau$spearman)))


enplusrf_rf_cau <- data.frame(spearman=filt_elnet_rf_cau_2_mets_0.1$rf_spearman, 
                              prediction=rep("EN+RF RF Only",length(filt_elnet_rf_cau_2_mets_0.1$rf_spearman)))

enplusrf_en_cau <- data.frame(spearman=filt_elnet_rf_cau_2_mets_0.1$spearman, 
                              prediction=rep("EN+RF EN Only",length(filt_elnet_rf_cau_2_mets_0.1$spearman)))


den_cau_2_mets <- rbind(rfonly_den_cau, enonly_den_cau, enplusrf_rf_cau, enplusrf_en_cau)

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("Random Forest", "Elastic NET", "Intersect"), no=c(424, 1439, 615)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (CAU 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine svr and knn here

#make density plot
#find the mean of each ML group
library(plyr)
mu <- ddply(den_cau_2_mets, "prediction", summarise, grp.mean=mean(spearman))

#Combine all density plot into one
# lwd = line thickness
ggplot(den_cau_2_mets, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Genes for RF and EN (CAU 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")


# make density plot for EN intersect RF HIS 2 METS
int_enrf_en_cau <- data.frame(spearman=filt_elnet_rf_cau_2_mets_0.1$spearman, 
                              prediction=rep("EN",length(filt_elnet_rf_cau_2_mets_0.1$spearman)))

int_enrf_rf_cau <- data.frame(spearman=filt_elnet_rf_cau_2_mets_0.1$rf_spearman, 
                              prediction=rep("RF",length(filt_elnet_rf_cau_2_mets_0.1$rf_spearman)))

den_int_enrf_cau <- rbind(int_enrf_en_cau, int_enrf_rf_cau)
#find the mean of each ML group
library(plyr)
mu <- ddply(den_int_enrf_cau, "prediction", summarise, grp.mean=mean(spearman))

# lwd = line thickness
ggplot(den_int_enrf_cau, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Intersect Genes for RF and EN (CAU 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")

#SVR
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
names(filt_svr_cau_2_mets) <- c("gene", "svr_spearman")

elnet_cau_2_mets_0.1 <- subset(elnet_cau_2_mets, spearman > 0.1)
filt_svr_cau_2_mets_0.1 <- subset(filt_svr_cau_2_mets, svr_spearman > 0.1)
filt_elnet_svr_cau_2_mets_0.1 <- inner_join(elnet_cau_2_mets_0.1, filt_svr_cau_2_mets_0.1, by = c("gene" = "gene"))
#plot, NOTE all models here are optimized
ggplot(filt_elnet_svr_cau_2_mets_0.1, aes(x=spearman, y=svr_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)") + 
  ylab("SVR") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_svr_cau_2_mets_0.1, x = "spearman", y = "svr_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "SVR", title = "Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#intersect before filtering off R > 0.1 for each of ML and EN
filt_elnet_svr_cau_2_mets <- inner_join(elnet_cau_2_mets, filt_svr_cau_2_mets, by = c("gene" = "gene"))
sub_elnet_svr_cau_2_mets <- subset(filt_elnet_svr_cau_2_mets, spearman > 0.1 | svr_spearman > 0.1) #then filter off where both are < 0.1

#Plot
ggplot(sub_elnet_svr_cau_2_mets, aes(x=spearman, y=svr_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)") + 
  ylab("SVR") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.4,1)) + ylim(c(-0.4,1)) + theme_bw()
#Then find the genes where EN or SVR is positive and the other negative 
#Do FUMA on top 10 of theses genes
svrpos_cau <- subset(sub_elnet_svr_cau_2_mets, spearman < 0)
enpos_s_cau <- subset(sub_elnet_svr_cau_2_mets, svr_spearman < 0)

##########
#get the genes that are not in the other. svr genes not in en. also en genes not in svr
#Do FUMA on top 10 of these genes
svronly_cau <- anti_join(filt_svr_cau_2_mets_0.1, elnet_cau_2_mets_0.1, by = c("gene" = "gene")) #Better!

enonly_s_cau <- anti_join(elnet_cau_2_mets_0.1, filt_svr_cau_2_mets_0.1,by = c("gene" = "gene"))

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("SVR", "Elastic NET", "Intersect"), no=c(410, 1411, 643)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (CAU 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here


#KNN
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
names(filt_knn_cau_2_mets) <- c("gene", "knn_spearman")

elnet_cau_2_mets_0.1 <- subset(elnet_cau_2_mets, spearman > 0.1)
filt_knn_cau_2_mets_0.1 <- subset(filt_knn_cau_2_mets, knn_spearman > 0.1)
filt_elnet_knn_cau_2_mets_0.1 <- inner_join(elnet_cau_2_mets_0.1, filt_knn_cau_2_mets_0.1, by = c("gene" = "gene"))
#plot Intersect where both ML R > 0.1, NOTE all models here are optimized
ggplot(filt_elnet_knn_cau_2_mets_0.1, aes(x=spearman, y=knn_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(0,1)) + ylim(c(0,1))

ggscatter(filt_elnet_knn_cau_2_mets_0.1, x = "spearman", y = "knn_spearman", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "KNN", title = "Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)",
          xlim = c(0, 1), ylim = c(0, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

#intersect before filtering off R > 0.1 for each of ML and EN
filt_elnet_knn_cau_2_mets <- inner_join(elnet_cau_2_mets, filt_knn_cau_2_mets, by = c("gene" = "gene"))
sub_elnet_knn_cau_2_mets <- subset(filt_elnet_knn_cau_2_mets, spearman > 0.1 | knn_spearman > 0.1) #then filter off where both are < 0.1

#Plot
ggplot(sub_elnet_knn_cau_2_mets, aes(x=spearman, y=knn_spearman)) + 
  ggtitle("Spearman Corr of Observed and Predicted Gene Expression (CAU to METS)") + 
  ylab("KNN") + xlab("Elastic Net") +
  geom_point(shape=1) + geom_abline(slope=1,intercept=0,col='blue') + xlim(c(-0.5,1)) + ylim(c(-0.5,1)) + theme_bw()
#Then find the genes where EN or KNN is positive and the other negative 
#Do FUMA on top 10 of theses genes
knnpos_cau <- subset(sub_elnet_knn_cau_2_mets, spearman < 0)
enpos_k_cau <- subset(sub_elnet_knn_cau_2_mets, knn_spearman < 0)

##########
#get the genes that are not in the other. svr genes not in en. also en genes not in svr
#Do FUMA on top 10 of these genes
knnonly_cau <- anti_join(filt_knn_cau_2_mets_0.1, elnet_cau_2_mets_0.1, by = c("gene" = "gene")) #Better!

enonly_k_cau <- anti_join(elnet_cau_2_mets_0.1, filt_knn_cau_2_mets_0.1,by = c("gene" = "gene"))

#find the name of those genes
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
###########
#Make bar plot with count of genes in only ML, EN, and EN intersect ML for genes R > 0.1
gcount <- data.frame(ml=c("KNN", "Elastic NET", "Intersect"), no=c(303, 1612, 442)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (CAU 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here

#Make bar plot with count of genes in only ML, EN+ML, and EN for genes R > 0.1
gcount <- data.frame(ml=c("EN-RF", "EN+RF", "RF-EN", "EN-SVR", "EN+SVR", "SVR-EN",
                          "EN-KNN", "EN+KNN", "KNN-EN"), no=c(1439,615,424,1411,643,410,1612,442,303)) #rfonly, enonly, 928
ggplot(data=gcount, aes(x=ml, y=no, fill=ml)) + geom_bar(stat="identity") + ylab("Gene Count") + xlab("Models") +
  ggtitle("Genes in Models (CAU 2 METS)") + geom_text(aes(label=no), vjust=1.6, color="white", size=3.5)+
  theme_minimal() # I will see if I can combine rf and knn here

#Make violin plot for all 3 population mesa to mets with EN
elnet_afa_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_AFA_2_METS.txt", header = T)
elnet_afa_2_mets$gene <- as.character(elnet_afa_2_mets$gene)
en_afa_2mets <- data.frame(spearman=elnet_afa_2_mets_0.1$spearman, prediction=rep("AFA_2_METS", length(elnet_afa_2_mets_0.1$spearman)))

elnet_his_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_HIS_2_METS.txt", header = T)
elnet_his_2_mets$gene <- as.character(elnet_his_2_mets$gene)
elnet_his_2_mets_0.1 <- subset(elnet_his_2_mets, spearman > 0.1)
en_his_2mets <- data.frame(spearman=elnet_his_2_mets_0.1$spearman, prediction=rep("HIS_2_METS", length(elnet_his_2_mets_0.1$spearman)))

elnet_cau_2_mets <- read.table(file = "Z:/data/mesa_models/spearman_CAU_2_METS.txt", header = T)
elnet_cau_2_mets$gene <- as.character(elnet_cau_2_mets$gene)
elnet_cau_2_mets_0.1 <- subset(elnet_cau_2_mets, spearman > 0.1)
en_cau_2mets <- data.frame(spearman=elnet_cau_2_mets_0.1$spearman, prediction=rep("CAU_2_METS", length(elnet_cau_2_mets_0.1$spearman)))

en_mesa_2_mets <- rbind(en_afa_2mets, en_his_2mets, en_cau_2mets)
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

library(Hmisc)
ggplot(en_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction, ylim(0,1))) + 
  geom_violin(trim = F) + geom_boxplot(width=0.2, color="black") + theme_minimal() + 
  ggtitle("Elastic Net (MESA to METS)") #+
#stat_summary(fun.y=mean, geom = "point", size=3, color="white")

#Make density plot for all 4 ML models with ML R > 0.1 CAU 2 METS
rf_0.1_cau <- data.frame(spearman=filt_rf_cau_2_mets_0.1$rf_spearman, 
                         prediction=rep("RF", length(filt_rf_cau_2_mets_0.1$rf_spearman)))

svr_0.1_cau <- data.frame(spearman=filt_svr_cau_2_mets_0.1$svr_spearman, 
                          prediction=rep("SVR", length(filt_svr_cau_2_mets_0.1$svr_spearman)))

knn_0.1_cau <- data.frame(spearman=filt_knn_cau_2_mets_0.1$knn_spearman, 
                          prediction=rep("KNN", length(filt_knn_cau_2_mets_0.1$knn_spearman)))

en_0.1_cau <- data.frame(spearman=elnet_cau_2_mets_0.1$spearman, 
                         prediction=rep("EN", length(elnet_cau_2_mets_0.1$spearman)))
den_4_ml_cau <- rbind(rf_0.1_cau, svr_0.1_cau, knn_0.1_cau, en_0.1_cau)


#find the mean of each ML group
library(plyr)
mu <- ddply(den_4_ml_cau, "prediction", summarise, grp.mean=mean(spearman))

# density plot for all ML in afa 2 mets
# lwd = line thickness
ggplot(den_4_ml_cau, aes(x = spearman, color=prediction, lwd=1.5)) + 
  scale_x_continuous(name = "Spearman Correlation > 0.1") + 
  ggtitle("Density plot of Genes for All Models (CAU 2 METS)") +
  geom_density()+ theme_classic() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=prediction),linetype="longdash", lwd=1) +
  geom_vline(xintercept = 0.5, size = 1, colour = "gold4", linetype = "dashed")











#Make violin plot for all 3 population mesa to mets with RF
rf_afa_2mets <- data.frame(spearman=filt_rf_afa_2_mets_0.1$rf_spearman, 
                           prediction=rep("AFA_2_METS", length(filt_rf_afa_2_mets_0.1$rf_spearman)))

rf_his_2mets <- data.frame(spearman=filt_rf_his_2_mets_0.1$rf_spearman, 
                           prediction=rep("HIS_2_METS", length(filt_rf_his_2_mets_0.1$rf_spearman)))

rf_cau_2mets <- data.frame(spearman=filt_rf_cau_2_mets_0.1$rf_spearman, 
                           prediction=rep("CAU_2_METS", length(filt_rf_cau_2_mets_0.1$rf_spearman)))

rf_mesa_2_mets <- rbind(rf_afa_2mets, rf_his_2mets, rf_cau_2mets)

ggplot(rf_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction)) + 
  geom_violin(trim=T) + geom_boxplot(width=0.2, color="black") +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.2), limits=c(0, 1.0))



#Make violin plot for all 3 population mesa to mets with SVR
svr_afa_2mets <- data.frame(spearman=filt_svr_afa_2_mets_0.1$svr_spearman, 
                           prediction=rep("AFA_2_METS", length(filt_svr_afa_2_mets_0.1$svr_spearman)))

svr_his_2mets <- data.frame(spearman=filt_svr_his_2_mets_0.1$svr_spearman, 
                           prediction=rep("HIS_2_METS", length(filt_svr_his_2_mets_0.1$svr_spearman)))

svr_cau_2mets <- data.frame(spearman=filt_svr_cau_2_mets_0.1$svr_spearman, 
                           prediction=rep("CAU_2_METS", length(filt_svr_cau_2_mets_0.1$svr_spearman)))

svr_mesa_2_mets <- rbind(svr_afa_2mets, svr_his_2mets, svr_cau_2mets)

ggplot(svr_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction, ylim(0,1))) + 
  geom_violin(trim=F) + geom_boxplot(width=0.2, color="black") + theme_minimal() + 
  ggtitle("SVR (MESA to METS)") + stat_summary(fun.y=mean, geom = "point", size=3, color="white")




#Make violin plot for all 3 population mesa to mets with KNN
knn_afa_2mets <- data.frame(spearman=filt_knn_afa_2_mets_0.1$knn_spearman, 
                            prediction=rep("AFA_2_METS", length(filt_knn_afa_2_mets_0.1$knn_spearman)))

knn_his_2mets <- data.frame(spearman=filt_knn_his_2_mets_0.1$knn_spearman, 
                            prediction=rep("HIS_2_METS", length(filt_knn_his_2_mets_0.1$knn_spearman)))

knn_cau_2mets <- data.frame(spearman=filt_knn_cau_2_mets_0.1$knn_spearman, 
                            prediction=rep("CAU_2_METS", length(filt_knn_cau_2_mets_0.1$knn_spearman)))

knn_mesa_2_mets <- rbind(knn_afa_2mets, knn_his_2mets, knn_cau_2mets)

ggplot(knn_mesa_2_mets, aes(x=prediction, y=spearman, color=prediction, fill=prediction, ylim(0,1))) + 
  geom_violin(trim=F) + geom_boxplot(width=0.2, color="black") + theme_minimal() + 
  ggtitle("KNN (MESA to METS)") + stat_summary(fun.y=mean, geom = "point", size=3, color="white")


#Gene Annotation. use it to find names of genes for FUMA
gv18 <- read.table(file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
gv18$gene_id <- as.character(gv18$gene_id)
gv18_intact <- gv18
for (i in 1:length(gv18$gene_id)){
  gv18$gene_id[i] <- gsub('\\.[0-9]+','',gv18$gene_id[i])
} #just to remove the decimal places in the gene_id

#Find genes where RF is positive and EN is negative across population
fuma_rf_pos <- rfpos[,c(1,3)]
fuma_rf_pos_his <- rfpos_his[,c(1,3)]
fuma_rf_pos_cau <-rfpos_cau[,c(1,3)]
#find there intersects
fuma_rf_afa_his <- inner_join(fuma_rf_pos, fuma_rf_pos_his, by = c("gene" = "gene")) # AFA and HIS
fuma_rf_afa_cau <- inner_join(fuma_rf_pos, fuma_rf_pos_cau, by = c("gene" = "gene")) # AFA and CAU
fuma_rf_his_cau <- inner_join(fuma_rf_pos_his, fuma_rf_pos_cau, by = c("gene" = "gene")) # HIS and CAU

#intersect all three pop
fuma_rf <- inner_join(fuma_rf_pos, fuma_rf_pos_his, by = c("gene" = "gene"))
fuma_rf <- inner_join(fuma_rf, fuma_rf_pos_cau, by = c("gene" = "gene"))


#Find genes with R > 0.1 found in RF and not in EN across population
fuma_rfonly <- rfonly
fuma_rfonly_his <- rfonly_his
fuma_rfonly_cau <- rfonly_cau
#find there intersects
fuma_rfonly_afa_his <- inner_join(fuma_rfonly, fuma_rfonly_his, by = c("gene"="gene")) # AFA and HIS
fuma_rfonly_afa_cau <- inner_join(fuma_rfonly, fuma_rfonly_cau, by = c("gene"="gene")) # AFA and CAU
fuma_rfonly_his_cau <- inner_join(fuma_rfonly_his, fuma_rfonly_cau, by = c("gene"="gene")) #HIS and CAU

#intersect for all three pop
fuma_only_rf <- inner_join(fuma_rfonly, fuma_rfonly_his, by = c("gene" = "gene"))
fuma_only_rf <- inner_join(fuma_only_rf, fuma_rfonly_cau, by = c("gene" = "gene")) #all three pop

fuma_gcode <- inner_join(fuma_only_rf, gv18, by = c("gene"="gene_id"))
