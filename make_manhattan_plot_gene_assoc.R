library(data.table)
library(dplyr)
#note that there is a bug in fread which makes it to skip some row

rfassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/rf_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt", header=T)
svrassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/svr_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt", header=T)
knnassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/knn_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt", header=T)


#make qqplot
#install.packages("qqman")
library(qqman)

qq(knnassoc$p, main="knn")
qq(svrassoc$p, main="svr")
qq(rfassoc$p, main="rf")

#qqnorm(knnassoc1$p) shows the pvalue came from a population that is normal distribution


#filter the pheno asscoiation result to have only ALL gene models with cv R2 > 0.01
#read in the best gird results for the ALL
all_rfgrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
all_svrgrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
all_knngrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)

all_rfgrid$Gene_ID <- as.character(all_rfgrid$Gene_ID)
all_svrgrid$Gene_ID <- as.character(all_svrgrid$Gene_ID)
all_knngrid$Gene_ID <- as.character(all_knngrid$Gene_ID)

#filter by cv R2 > 0.01
all_rfgrid <- subset(all_rfgrid, CV_R2 > 0.01)
all_rfgrid <- all_rfgrid[,c(1,2,3)]# keep only gene_id, gene_name, and cv_r2
all_svrgrid <- subset(all_svrgrid, CV_R2 > 0.01)
all_svrgrid <- all_svrgrid[,c(1,2,3)]
all_knngrid <- subset(all_knngrid, CV_R2 > 0.01)
all_knngrid <- all_knngrid[,c(1,2,3)]

#retain only genes with CV R2 > 0.01 in pheno assoc
rfassoc$gene <- as.character(rfassoc$gene)
svrassoc$gene <- as.character(svrassoc$gene)
knnassoc$gene <- as.character(knnassoc$gene)

rfassoc_cvr0.01 <- inner_join(all_rfgrid, rfassoc, by = c("Gene_ID" = "gene"))
svrassoc_cvr0.01 <- inner_join(all_svrgrid, svrassoc, by = c("Gene_ID" = "gene"))
knnassoc_cvr0.01 <- inner_join(all_knngrid, knnassoc, by = c("Gene_ID" = "gene"))

#ggplot after filtering
qq(rfassoc_cvr0.01$p, main="rf filtered by cv r2")
qq(svrassoc_cvr0.01$p, main="svr filtered by cv r2")
qq(knnassoc_cvr0.01$p, main="knn filtered by cv r2")
rf_pheno_qq_filtered


#Make manhattan plot for the gene pheno assoc results
library(qqman)
#see example table for manhattan plot
eggwas <- gwasResults #this comes with qqman

#merge the gene pheno assoc results with gencode annotation
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode1 <- gencode[,c(1,2,3,4)]

rfassocgen <- inner_join(gencode1, rfassoc, by = c("gene_id" = "gene"))
rfassocgen <- rfassocgen[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svrassocgen <- inner_join(gencode1, svrassoc, by = c("gene_id" = "gene"))
svrassocgen <- svrassocgen[,c(1,3,4,7)] #drop gene_id

knnassocgen <- inner_join(gencode1, knnassoc, by = c("gene_id" = "gene"))
knnassocgen <- knnassocgen[,c(1,3,4,7)] #drop gene_id

#take each chrom out, sort it by "start", and rbind them
rfassocman <- NULL
svrassocman <- NULL
knnassocman <- NULL

#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rfassocgen, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="rf with cv_R2 filtering", annotatePval=0.0001,
          col=c("red", "blue"))

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svrassocgen, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="svr without cv_R2 filtering", annotatePval=0.001,
          col=c("red", "blue"))

#knn
for (i in 1:22){
  i <- as.character(i)
  a <- subset(knnassocgen, chr==i)
  a <- a[order(a$start),]
  knnassocman <- rbind(knnassocman,a)
}
names(knnassocman) <- c("CHR", "SNP", "BP", "P")
knnassocman$CHR <- as.numeric(knnassocman$CHR)
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="knn without cv_R2 filtering", 
          annotatePval=0.0001, col=c("red", "blue"))


