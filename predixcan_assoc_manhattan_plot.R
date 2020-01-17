library(data.table)
library(dplyr)
#note that there is a bug in fread which makes it to skip some row for huge files

enassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/en_cau_rankplt5_association.txt", header=T)
enassoc$gene <- as.character(enassoc$gene)
for (i in 1:nrow(enassoc)){
  enassoc$gene[i] <- gsub('\\.[0-9]+','',enassoc$gene[i])
} #just to remove the decimal places in the gene_id


#make qqplot
#install.packages("qqman")
library(qqman)

qq(enassoc$p, main="en")

#qqnorm(knnassoc1$p) shows the pvalue came from a population that is normal distribution


#filter the pheno asscoiation result to have only ALL gene models with cv R2 > 0.01
#read in the best en results for the ALL
#file does not have header, so read in ones with header and copy the header
en_cau <- fread(file="Z:/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
en_all <- fread(file="Z:/no_header_MESA.ALL.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=F)
header <- colnames(en_cau)
names(en_all) <- header

en_all$gene_id <- as.character(en_all$gene_id)
for (i in 1:nrow(en_all)){
  en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
} #just to remove the decimal places in the gene_id


en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)

#filter by cv R2 > 0.01
en_all <- subset(en_all, cv_R2_avg > 0.01)
en_all <- en_all[,c(1,2,10)]# keep only gene_id, gene_name, and cv_r2


#retain only genes with CV R2 > 0.01 in pheno assoc

enassoc_cvr0.01 <- inner_join(en_all, enassoc, by = c("gene_id" = "gene"))

#write out the filtered assoc
fwrite(enassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/en_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")

#ggplot after filtering
qq(enassoc_cvr0.01$p, main="en filtered by cv r2")


#Make manhattan plot for the gene pheno assoc results
library(qqman)
#see example table for manhattan plot
eggwas <- gwasResults #this comes with qqman

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

enassocgen <- inner_join(gencode1, enassoc_cvr0.01, by = c("gene_id" = "gene_id"))
enassocgen <- enassocgen[,c(1,3,4,9)] #drop gene_id and use gene_name to help annotation in the plot

#take each chrom out, sort it by "start", and rbind them
enassocman <- NULL

#en
for (i in 1:22){
  i <- as.character(i)
  a <- subset(enassocgen, chr==i)
  a <- a[order(a$start),]
  enassocman <- rbind(enassocman,a)
}
names(enassocman) <- c("CHR", "SNP", "BP", "P")
enassocman$CHR <- as.numeric(enassocman$CHR)
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="predixcan with cv_R2 filtering corrected", annotatePval=0.0005,
          col=c("red", "blue"))

