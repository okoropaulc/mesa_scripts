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
