library(data.table)
library(dplyr)
#note that there is a bug in fread which makes it to skip some row for huge files

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
all_rfgrid$CV_R2 <- as.numeric(all_rfgrid$CV_R2)
all_svrgrid$Gene_ID <- as.character(all_svrgrid$Gene_ID)
all_svrgrid$CV_R2 <- as.character(all_svrgrid$CV_R2)
all_knngrid$Gene_ID <- as.character(all_knngrid$Gene_ID)
all_knngrid$CV_R2 <- as.numeric(all_knngrid$CV_R2)

#first remove decimals from the gene id. the number of genes=9623 is same for all algs
for (i in 1:9623){
  all_rfgrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_rfgrid$Gene_ID[i])
  all_svrgrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_svrgrid$Gene_ID[i])
  all_knngrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_knngrid$Gene_ID[i])
} #just to remove the decimal places in the gene_id

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

#first remove decimals from the gene id. the number of genes=9623 is same for all algs
for (i in 1:9623){
  rfassoc$gene[i] <- gsub('\\.[0-9]+','',rfassoc$gene[i])
  svrassoc$gene[i] <- gsub('\\.[0-9]+','',svrassoc$gene[i])
  knnassoc$gene[i] <- gsub('\\.[0-9]+','',knnassoc$gene[i])
} #just to remove the decimal places in the gene_id

rfassoc_cvr0.01 <- inner_join(all_rfgrid, rfassoc, by = c("Gene_ID" = "gene"))
svrassoc_cvr0.01 <- inner_join(all_svrgrid, svrassoc, by = c("Gene_ID" = "gene"))
knnassoc_cvr0.01 <- inner_join(all_knngrid, knnassoc, by = c("Gene_ID" = "gene"))

#write out the filtered assoc
fwrite(rfassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/rf_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/svr_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/knn_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")


#ggplot after filtering
qq(rfassoc_cvr0.01$p, main="rf filtered by cv r2")
qq(svrassoc_cvr0.01$p, main="svr filtered by cv r2")
qq(knnassoc_cvr0.01$p, main="knn filtered by cv r2")
#rf_pheno_qq_filtered


#Make manhattan plot for the gene pheno assoc results
library(qqman)
#see example table for manhattan plot
#eggwas <- gwasResults #this comes with qqman

#manhattan for assoc filtered by cv r2
#merge the gene pheno assoc results with gencode annotation
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

#gencode1 <- subset(gencode, gene_type=="protein_coding")
#gencode1$chr <- as.numeric(gencode1$chr)
chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
gencode1 <- gencode[,c(1,2,3,4)]

rfassocgen <- inner_join(gencode1, rfassoc_cvr0.01, by = c("gene_id" = "Gene_ID"))
rfassocgen <- rfassocgen[,c(1,3,4,9)] #drop gene_id and use gene_name to help annotation in the plot

svrassocgen <- inner_join(gencode1, svrassoc_cvr0.01, by = c("gene_id" = "Gene_ID"))
svrassocgen <- svrassocgen[,c(1,3,4,9)] #drop gene_id

knnassocgen <- inner_join(gencode1, knnassoc_cvr0.01, by = c("gene_id" = "Gene_ID"))
knnassocgen <- knnassocgen[,c(1,3,4,9)] #drop gene_id

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
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="rf with cv_R2 filtering corrected", annotatePval=0.0005,
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
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="svr with cv_R2 filtering corrected", annotatePval=0.0005,
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
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="knn with cv_R2 filtering corrected", 
          annotatePval=0.0005, col=c("red", "blue"))









# #merge the gene pheno assoc results with gencode annotation
# gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
# gencode$gene_id <- as.character(gencode$gene_id)
# gencode1 <- subset(gencode, gene_type=="protein_coding")
# #gencode1$chr <- as.numeric(gencode1$chr)
# chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
# #gencode1 <- gencode[,c(1,2,3,4)]
# 
# rfassocgen <- inner_join(gencode1, rfassoc, by = c("gene_id" = "gene"))
# rfassocgen <- rfassocgen[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot
# 
# svrassocgen <- inner_join(gencode1, svrassoc, by = c("gene_id" = "gene"))
# svrassocgen <- svrassocgen[,c(1,3,4,7)] #drop gene_id
# 
# knnassocgen <- inner_join(gencode1, knnassoc, by = c("gene_id" = "gene"))
# knnassocgen <- knnassocgen[,c(1,3,4,7)] #drop gene_id
# 
# #take each chrom out, sort it by "start", and rbind them
# rfassocman <- NULL
# svrassocman <- NULL
# knnassocman <- NULL
# 
# #rf
# for (i in 1:22){
#   i <- as.character(i)
#   a <- subset(rfassocgen, chr==i)
#   a <- a[order(a$start),]
#   rfassocman <- rbind(rfassocman,a)
# }
# names(rfassocman) <- c("CHR", "SNP", "BP", "P")
# rfassocman$CHR <- as.numeric(rfassocman$CHR)
# manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="rf with cv_R2 filtering", annotatePval=0.0001,
#           col=c("red", "blue"))
# 
# #svr
# for (i in 1:22){
#   i <- as.character(i)
#   a <- subset(svrassocgen, chr==i)
#   a <- a[order(a$start),]
#   svrassocman <- rbind(svrassocman,a)
# }
# names(svrassocman) <- c("CHR", "SNP", "BP", "P")
# svrassocman$CHR <- as.numeric(svrassocman$CHR)
# manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="svr without cv_R2 filtering", annotatePval=0.001,
#           col=c("red", "blue"))
# 
# #knn
# for (i in 1:22){
#   i <- as.character(i)
#   a <- subset(knnassocgen, chr==i)
#   a <- a[order(a$start),]
#   knnassocman <- rbind(knnassocman,a)
# }
# names(knnassocman) <- c("CHR", "SNP", "BP", "P")
# knnassocman$CHR <- as.numeric(knnassocman$CHR)
# manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="knn without cv_R2 filtering", 
#           annotatePval=0.0001, col=c("red", "blue"))
# 
# 
