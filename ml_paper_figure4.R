#plot figure 4 in ML Paper
library(data.table)
library(ggplot2)
library(dplyr)

#read in the assoc again so R can treat the pvalue as numeric
assoc_en <- fread(file="Z:/data/twas_mesa/en_rank_hdl_assoc.txt",header=T)
assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rank_hdl_assoc.txt",header=T)
assoc_svr <- fread(file="Z:/data/twas_mesa/svr_rank_hdl_assoc.txt",header=T)
assoc_knn <- fread(file="Z:/data/twas_mesa/knn_rank_hdl_assoc.txt",header=T)


library(qqman)
p1 <- qq(assoc_en$p, main="Elastic Net Rank_HDL", ylim=c(0,14), cex=2, las=1)
p2 <- qq(assoc_rf$p, main="Random Forest Rank_HDL", ylim=c(0,14), cex=2, las=1)
p3 <- qq(assoc_svr$p, main="Support Vector Rank_HDL", ylim=c(0,14), cex=2, las=1)
p4 <- qq(assoc_knn$p, main="KNN Rank_HDL", ylim=c(0,14), cex=2, las=1)

en <- data.frame(p=assoc_en$p)
en <- as.data.frame(-log10(en$p))
names(en) <- "p"
ggplot(en, aes(sample=p)) + stat_qq() + stat_qq_line()

#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2)

#merge the gene pheno assoc results with gencode annotation
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

gencode1 <- gencode[,c(1,2,3,4)]

en_man <- inner_join(gencode1, assoc_en, by = c("gene_id" = "gene"))
en_man <- en_man[,c(1,3,4,7)]#drop gene_id and use gene_name to help annotation in the plot

rf_man <- inner_join(gencode1, assoc_rf, by = c("gene_id" = "gene"))
rf_man <- rf_man[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svr_man <- inner_join(gencode1, assoc_svr, by = c("gene_id" = "gene"))
svr_man <- svr_man[,c(1,3,4,7)] #drop gene_id

knn_man <- inner_join(gencode1, assoc_knn, by = c("gene_id" = "gene"))
knn_man <- knn_man[,c(1,3,4,7)] #drop gene_id


#take each chrom out, sort it by "start", and rbind them
enassocman <- NULL
rfassocman <- NULL
svrassocman <- NULL
knnassocman <- NULL

#en
for (i in 1:22){
  i <- as.character(i)
  a <- subset(en_man, chr==i)
  a <- a[order(a$start),]
  enassocman <- rbind(enassocman,a)
}
names(enassocman) <- c("CHR", "SNP", "BP", "P")
enassocman$CHR <- as.numeric(enassocman$CHR)
p5 <- manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net HDL", suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=2, cex.axis=1.5, cex.lab=1)

#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
p6 <- manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest HDL", suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=2, cex.axis=1.5, cex.lab=1)

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
p7 <- manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector HDL", suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=2, cex.axis=1.5, cex.lab=1)

#knn
for (i in 1:22){
  i <- as.character(i)
  a <- subset(knn_man, chr==i)
  a <- a[order(a$start),]
  knnassocman <- rbind(knnassocman,a)
}
names(knnassocman) <- c("CHR", "SNP", "BP", "P")
knnassocman$CHR <- as.numeric(knnassocman$CHR)
knnassocman$P <- as.numeric(knnassocman$P)
p8 <- manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN HDL", suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=2, cex.axis=1.5, cex.lab=1)

grid.arrange(p5,p6,p7,p8,nrow=2)
