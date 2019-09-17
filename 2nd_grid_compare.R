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

ggscatter(elrf, x = "elnet", y = "rf", add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "RF Optimized", title = "Spearman Corr of Observed and Predicted Gene Expression (AFA to METS)",
          xlim = c(-1, 1), ylim = c(-1, 1)) + geom_abline(intercept = 0, slope = 1, color="blue")

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
 