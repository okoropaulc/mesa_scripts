library(dplyr)
library(ggplot2)
"%&%" <- function(a,b) paste(a,b, sep = "")


gene_annotation_file <- read.table(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
#gene_annotation_file$gene_type <- as.character(gene_annotation_file$gene_type)
gene_annotation_file <- subset(gene_annotation_file, gene_type == "protein_coding")
gene_annotation_file$chr <- as.character(gene_annotation_file$chr)
chr22 <- subset(gene_annotation_file, chr=="22")
chr22$gene_id <- as.character(chr22$gene_id)
chr22 <- data.frame(gene_id=chr22[,2])

# Max_evals=30
chr1chk1 <- read.table(file="Z:/data/paper_hyperopt/afa_en_hyperopt_chr1_chunk1.txt", header=T, sep="\t")
chr1chk2 <- read.table(file="Z:/data/paper_hyperopt/afa_en_hyperopt_chr1_chunk2.txt", header=T, sep="\t")
chr1chk3 <- read.table(file="Z:/data/paper_hyperopt/afa_en_hyperopt_chr1_chunk3.txt", header=T, sep="\t")
chr1chk4 <- read.table(file="Z:/data/paper_hyperopt/afa_en_hyperopt_chr1_chunk4.txt", header=T, sep="\t")
chr1chk5 <- read.table(file="Z:/data/paper_hyperopt/afa_en_hyperopt_chr1_chunk5.txt", header=T, sep="\t")

chr1 <- rbind(chr1chk1, chr1chk2, chr1chk3, chr1chk4, chr1chk5)
chr1$X <- NULL
pen <- chr1[,c(1,34)]


en_afa <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header=T)
ren <- en_afa[,c(1,10)]


library(dplyr)

data <- inner_join(pen, ren, by = c("gene_id"="gene_id"))
names(data) <- c("gene", "P", "R")

#remove outliers
data <- subset(data, (P>-1 & R>-1))

library(ggplot2)

ggplot(data = data, aes(x=R, y=P)) + geom_point()


#KNN
chr1chk1 <- read.table(file="Z:/data/paper_hyperopt/afa_knn_hyperopt_chr1_chunk1.txt", header=T, sep="\t")
chr1chk2 <- read.table(file="Z:/data/paper_hyperopt/afa_knn_hyperopt_chr1_chunk2.txt", header=T, sep="\t")
chr1chk3 <- read.table(file="Z:/data/paper_hyperopt/afa_knn_hyperopt_chr1_chunk3.txt", header=T, sep="\t")
chr1chk4 <- read.table(file="Z:/data/paper_hyperopt/afa_knn_hyperopt_chr1_chunk4.txt", header=T, sep="\t")
chr1chk5 <- read.table(file="Z:/data/paper_hyperopt/afa_knn_hyperopt_chr1_chunk5.txt", header=T, sep="\t")

chr1 <- rbind(chr1chk1, chr1chk2, chr1chk3, chr1chk4, chr1chk5)
chr1$X <- NULL
hknn <- chr1[,c(1,34)]



knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_chr1_full.txt", header = T)
gknn <- knn_afa[,c(1,3)]

# 
# data <- inner_join(hknn, gknn, by = c("gene_id"="Gene_ID"))
# names(data) <- c("gene", "H", "G")
# 
# ggplot(data = data, aes(x=G, y=H)) + geom_point()
# 
# 
# 
# #EN R Python Compare
# pchr22 <- read.table(file="Z:/data/mesa_models/en_R_Python_compare/AFA_en_py_chr22.txt", header=T)
# #rchr22 <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header=T)
# 
# pchr22 <- pchr22[,c(1,4)]
# rchr22 <- en_afa[,c(1,10)]
# 
# data <- inner_join(pchr22, rchr22, by = c("gene_id"="gene_id"))
# names(data) <- c("gene", "P", "R")
# 
# ggplot(data = data, aes(x=R, y=P)) + geom_point()



#merge all the chunks
"%&%" <- function(a,b) paste(a,b, sep = "")

#EN merge

pchr <- NULL

for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    pchr <- rbind(pchr, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
pchr$X <- NULL
pen <- pchr[,c(1,34)]
en_afa <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header=T)
ren <- en_afa[,c(1,10)]

library(dplyr)

data <- inner_join(pen, ren, by = c("gene_id"="gene_id"))
names(data) <- c("gene", "H", "G")
#remove outliers
data <- subset(data, (H>-1 & G>-1)) #243
library(ggplot2)
ggplot(data = data, aes(x=G, y=H)) + geom_point() + theme_bw(24) +ggtitle("EN") + ylim(-1,1) + xlim(-0.5,1)



#SVR
svr <- NULL

for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    svr <- rbind(svr, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}

svr$X <- NULL
hsvr <- svr[,c(1,34)]

svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header = T)
gsvr <- svr_afa[,c(1,3)]
data <- inner_join(hsvr, gsvr, by = c("gene_id"="Gene_ID"))
names(data) <- c("gene", "H", "G")
data <- subset(data, (H>-1 & G>-1)) #230
ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24) +ggtitle("SVR") + ylim(-1,1) + xlim(-0.5,1)


#KNN
knn <- NULL

for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    knn <- rbind(knn, read.table(file="Z:/data/paper_hyperopt/max_evals_30/afa_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}

knn$X <- NULL
hknn <- knn[,c(1,34)]

knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header = T)
gknn <- knn_afa[,c(1,3)]
data <- inner_join(hknn, gknn, by = c("gene_id"="Gene_ID"))
names(data) <- c("gene", "H", "G")
data <- subset(data, (H>-1 & G>-1)) #243
ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24) +ggtitle("KNN") + ylim(-1,1) + xlim(-0.5,1)


#compare the hyperopts of EN against SVR or KNN 
#SVR
ensvr <- inner_join(pen, hsvr, by = c("gene_id"="gene_id"))
names(ensvr) <- c("gene", "EN", "SVR")
ensvr <- subset(ensvr, (EN>-0.5 & SVR>-0.5)) #9441
ggplot(data = ensvr, aes(x=EN, y=SVR)) + geom_point() + geom_point() + theme_bw(24) + 
  geom_abline(intercept=0, slope=1, color="blue") + 
  geom_smooth(method="lm", color="red", lwd=0.5)

#KNN
enknn <- inner_join(pen, hknn, by = c("gene_id"="gene_id"))
names(enknn) <- c("gene", "EN", "KNN")
enknn <- subset(enknn, (EN>-1 & KNN>-1)) #9620
ggplot(data = enknn, aes(x=EN, y=KNN)) + geom_point() + geom_point() + theme_bw(24)




# ##  RF, chr22
# rf <- NULL
# 
# for (chrom in 22:22) {
#   no <- as.character(chrom)
#   for (chunk in 1:22) {
#     rf <- rbind(rf, read.table(file="Z:/data/paper_hyperopt/RF/afa_rf_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
#   }
# }
# 
# rf$X <- NULL
# hrf <- rf[,c(1,34)]
# 
# rf_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_rf_all_chrom.txt", header = T)
# grf <- rf_afa[,c(1,3)]
# 
# 
# data <- inner_join(hrf, grf, by = c("gene_id"="Gene_ID"))
# names(data) <- c("gene", "H", "G")
# 
# data <- subset(data, (H>-1 & G>-1)) #243
# ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24)
# 
# #compare the hyperopts of EN against RF
# enrf <- inner_join(pen, hrf, by = c("gene_id"="gene_id"))
# names(enrf) <- c("gene", "EN", "RF")
# enrf <- subset(enrf, (EN>-1 & RF>-1)) #9620
# ggplot(data = enrf, aes(x=EN, y=RF)) + geom_point() + geom_point() + theme_bw(24)



#Max_evals = 100

#EN merge
pchr <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    pchr <- rbind(pchr, read.table(file="Z:/data/paper_hyperopt/max_evals_100/afa_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
pchr$X <- NULL
pen <- pchr[,c(1,104)]

ren <- en_afa[,c(1,10)]

data <- inner_join(pen, ren, by = c("gene_id"="gene_id"))
names(data) <- c("gene", "H", "G")
#remove outliers
data <- subset(data, (H>-1 & G>-1)) #243
library(ggplot2)
ggplot(data = data, aes(x=G, y=H)) + geom_point() + theme_bw(24) +ggtitle("EN") + ylim(-1,1) + xlim(-0.5,1)


#SVR
svr <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    svr <- rbind(svr, read.table(file="Z:/data/paper_hyperopt/max_evals_100/afa_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}

svr$X <- NULL
hsvr <- svr[,c(1,104)]

svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header = T)
gsvr <- svr_afa[,c(1,3)]

data <- inner_join(hsvr, gsvr, by = c("gene_id"="Gene_ID"))
names(data) <- c("gene", "H", "G")
data <- subset(data, (H>-1 & G>-1)) #230
ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24) +ggtitle("SVR") + ylim(-1,1) + xlim(-0.5,1)


#KNN
knn <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:5) {
    knn <- rbind(knn, read.table(file="Z:/data/paper_hyperopt/max_evals_100/afa_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}

knn$X <- NULL
hknn <- knn[,c(1,104)]

knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header = T)
gknn <- knn_afa[,c(1,3)]
data <- inner_join(hknn, gknn, by = c("gene_id"="Gene_ID"))
names(data) <- c("gene", "H", "G")
data <- subset(data, (H>-1 & G>-1)) #243
ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24) +ggtitle("KNN") + ylim(-1,1) + xlim(-0.5,1)


#compare the hyperopts of EN against SVR or KNN 
#SVR
ensvr <- inner_join(pen, hsvr, by = c("gene_id"="gene_id"))
names(ensvr) <- c("gene", "EN", "SVR")
ensvr <- subset(ensvr, (EN>-1 & SVR>-1)) #9443
ggplot(data = ensvr, aes(x=EN, y=SVR)) + geom_point() + theme_bw(24)

#KNN
enknn <- inner_join(pen, hknn, by = c("gene_id"="gene_id"))
names(enknn) <- c("gene", "EN", "KNN")
enknn <- subset(enknn, (EN>-1 & KNN>-1)) #9616
ggplot(data = enknn, aes(x=EN, y=KNN)) + geom_point() + theme_bw(24)



# ###chr22 Chunk1-22
# 
# #EN
# en_afa <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header=T)
# ren <- en_afa[,c(1,10)]
# 
# pchr <- NULL
# for (chrom in 22:22) {
#   no <- as.character(chrom)
#   for (chunk in 1:22) {
#     pchr <- rbind(pchr, read.table(file="Z:/data/paper_hyperopt/EN/afa_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
#   }
# }
# pchr$X <- NULL
# pen <- pchr[,c(1,1004)]
# 
# ren <- en_afa[,c(1,10)]
# 
# library(dplyr)
# 
# data <- inner_join(pen, ren, by = c("gene_id"="gene_id"))
# names(data) <- c("gene", "P", "R")
# 
# #remove outliers
# data <- subset(data, (P>-1 & R>-1)) #9602
# library(ggplot2)
# 
# ggplot(data = data, aes(x=R, y=P)) + geom_point() + theme_bw(24)
# 
# 
# #SVR
# svr <- NULL
# 
# for (chrom in 22:22) {
#   no <- as.character(chrom)
#   for (chunk in 1:22) {
#     svr <- rbind(svr, read.table(file="Z:/data/paper_hyperopt/SVR/afa_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
#   }
# }
# 
# svr$X <- NULL
# hsvr <- svr[,c(1,104)]
# 
# svr_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header = T)
# gsvr <- svr_afa[,c(1,3)]
# 
# 
# data <- inner_join(hsvr, gsvr, by = c("gene_id"="Gene_ID"))
# names(data) <- c("gene", "H", "G")
# 
# data <- subset(data, (H>-1 & G>-1)) #230
# ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24)
# 
# 
# #KNN
# 
# knn <- NULL
# 
# for (chrom in 1:22) {
#   no <- as.character(chrom)
#   for (chunk in 1:5) {
#     knn <- rbind(knn, read.table(file="Z:/data/paper_hyperopt/afa_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
#   }
# }
# 
# knn$X <- NULL
# hknn <- knn[,c(1,34)]
# 
# knn_afa <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header = T)
# gknn <- knn_afa[,c(1,3)]
# 
# 
# data <- inner_join(hknn, gknn, by = c("gene_id"="Gene_ID"))
# names(data) <- c("gene", "H", "G")
# 
# data <- subset(data, (H>-1 & G>-1)) #243
# ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24)




#Max evals = 1000

#EN
en <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:22) {
    en <- rbind(en, read.table(file="Z:/data/paper_hyperopt/EN/afa_en_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
en$X <- NULL
en <- en[,c(1,1004)]

gen <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header=T)
gen <- gen[,c(1,10)]
data <- inner_join(en, gen, by = c("gene_id"="gene_id"))
names(data) <- c("gene", "H", "G")
data <- subset(data, (H>-1 & G>-1)) #233
ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24) + ggtitle("EN") + ylim(-1,1) + xlim(-0.5,1)

#SVR
svr <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:22) {
    svr <- rbind(svr, read.table(file="Z:/data/paper_hyperopt/SVR/afa_svr_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
svr$X <- NULL
svr <- svr[,c(1,1004)]

gsvr <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_svr_all_chrom.txt", header = T)
gsvr <- gsvr[,c(1,3)]
data <- inner_join(svr, gsvr, by = c("gene_id"="Gene_ID"))
names(data) <- c("gene", "H", "G")
data <- subset(data, (H>-1 & G>-1)) #232
ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24) +ggtitle("SVR") + ylim(-1,1) + xlim(-0.5,1)


#KNN
knn <- NULL
for (chrom in 22:22) {
  no <- as.character(chrom)
  for (chunk in 1:22) {
    knn <- rbind(knn, read.table(file="Z:/data/paper_hyperopt/KNN/afa_knn_hyperopt_chr" %&% no %&% "_chunk" %&% chunk %&% ".txt", header=T, sep="\t"))
  }
}
knn$X <- NULL
knn <- knn[,c(1,1004)]

gknn <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_knn_all_chrom.txt", header = T)
gknn <- gknn[,c(1,3)]
data <- inner_join(knn, gknn, by = c("gene_id"="Gene_ID"))
names(data) <- c("gene", "H", "G")
data <- subset(data, (H>-1 & G>-1)) #243
ggplot(data = data, aes(x=G, y=H)) + geom_point() + geom_point() + theme_bw(24) +ggtitle("KNN") + ylim(-1,1) + xlim(-0.5,1)

#SVR
ensvr <- inner_join(en, svr, by = c("gene_id"="gene_id"))
names(ensvr) <- c("gene", "EN", "SVR")
ensvr <- subset(ensvr, (EN>-0.5 & SVR>-0.5)) #9441
ggplot(data = ensvr, aes(x=EN, y=SVR)) + geom_point() + geom_point() + theme_bw(24) + 
  geom_abline(intercept=0, slope=1, color="blue") + 
  geom_smooth(method="lm", color="red", lwd=0.5)


#KNN
enknn <- inner_join(en, knn, by = c("gene_id"="gene_id"))
names(enknn) <- c("gene", "EN", "KNN")
enknn <- subset(enknn, (EN>-0.5 & KNN>-0.5)) #9441
ggplot(data = enknn, aes(x=EN, y=KNN)) + geom_point() + geom_point() + theme_bw(24) + 
  geom_abline(intercept=0, slope=1, color="blue") + 
  geom_smooth(method="lm", color="red", lwd=0.5)
