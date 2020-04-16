library(dplyr)
library(ggplot2)

###########               ALL

en_all <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", header=T)#it has no header
en_all <- en_all[,c(1,2,10)]
library(tidyverse)
en_all <- drop_na(en_all) #remove NA
en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)
en_all <- subset(en_all, cv_R2_avg > 0.01)
en_all$gene_id <- as.character(en_all$gene_id)
names(en_all) <- c("gene_id","gene_name","cvr2")


rf_all <- read.table(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
rf_all$CV_R2 <- as.numeric(rf_all$CV_R2)
rf_all <- subset(rf_all, CV_R2 > 0.01)
rf_all$Gene_ID <- as.character(rf_all$Gene_ID)
rf_all <- rf_all[,c(1,2,3)]
names(rf_all) <- c("gene_id","gene_name","cvr2")

en_only <- anti_join(en_all, rf_all, by = c("gene_id"="gene_id"))
rf_only <- anti_join(rf_all, en_all, by = c("gene_id"="gene_id"))

en_rf <- inner_join(en_all, rf_all, by = c("gene_id"="gene_id"))

ggplot(en_rf, aes(x=cvr2.x, y=cvr2.y)) + geom_point()


################## Check the comparison on the TWAS of ALL

library(data.table)
library(dplyr)
"%&%" <- function(a,b) paste(a,b, sep = "")
library(ggplot2)

#gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)

#cetp <- "ENSG00000087237"
#st8sia4 <- "ENSG00000113532"

pop <- "" #empty string indicates ALL
trait <- "hdl"

#for (trait in c("hdl", "ldl", "chol", "trig", "plt5")){}

#for (pop in c("AFA", "CAU", "HIS")){}

en <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "en_rank_" %&% trait %&% "_assoc.txt", header=T)
rf <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "rf_rank_" %&% trait %&% "_assoc.txt", header=T)
svr <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "svr_rank_" %&% trait %&% "_assoc.txt", header=T)
knn <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "knn_rank_" %&% trait %&% "_assoc.txt", header=T)



gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id


gencode1 <- gencode[,c(1,2,3,4)]

en_man <- inner_join(gencode1, en, by = c("gene_id" = "gene"))
en_man <- en_man[,c(1,3,4,7)]#drop gene_id and use gene_name to help annotation in the plot

rf_man <- inner_join(gencode1, rf, by = c("gene_id" = "gene"))
rf_man <- rf_man[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svr_man <- inner_join(gencode1, svr, by = c("gene_id" = "gene"))
svr_man <- svr_man[,c(1,3,4,7)] #drop gene_id

knn_man <- inner_join(gencode1, knn, by = c("gene_id" = "gene"))
knn_man <- knn_man[,c(1,3,4,7)] #drop gene_id


####check unique twas genes
en_twas_only <- anti_join(en_man, rf_man, by = c("gene_name"="gene_name"))
rf__twas_only <- anti_join(rf_man, en_man, by = c("gene_name"="gene_name"))
