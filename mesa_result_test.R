library(dplyr)

old_afa_chr14 <- read.table(file = "Z:/mesa_models/completed_old_train_results/new_AFA_model_nested_cv_chr14_model_summaries_10_peer_3pcs.txt", header = TRUE)
elnet_col <- colnames(old_afa_chr14)
afa_chr14$gene_id <- as.character(afa_chr14$gene_id)
afa_chr14 <- subset(afa_chr14,afa_chr14$gene_type == "protein_coding")

afa_chr13 <- read.table(file = "/home/paul/mesa_models/completed_old_train_results/new_AFA_model_nested_cv_chr13_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr13$gene_id <- as.character(afa_chr13$gene_id)
afa_chr13 <- subset(afa_chr13,afa_chr13$gene_type == "protein_coding")

prot_afa_chr13 <- read.table(file = "Z:/mesa_models/AFA.chr13.txt", header = FALSE)
colnames(prot_afa_chr13) <- elnet_col
prot_afa_chr13[,25] <- NULL
prot_afa_chr13 <- subset(prot_afa_chr13, prot_afa_chr13$cv_R2_avg >= 0.01)
prot_afa_chr13$gene_name <- as.character(prot_afa_chr13$gene_name)

prot_afa_chr14 <- read.table(file = "Z:/mesa_models/AFA.chr14.txt", header = FALSE)
colnames(prot_afa_chr14) <- elnet_col
prot_afa_chr14[,25] <- NULL
prot_afa_chr14 <- subset(prot_afa_chr14, prot_afa_chr14$cv_R2_avg >= 0.01)
prot_afa_chr14$gene_name <- as.character(prot_afa_chr14$gene_name)

AFA_orginal_R2 <- read.table(file = "Z:/mesa_models/AFA_model_summary_stats.txt", header = TRUE)
colnames(AFA_orginal_R2) <- elnet_col
AFA_orginal_R2 <- subset(AFA_orginal_R2,AFA_orginal_R2$gene_type == "protein_coding")
AFA_orginal_R2[,25] <- NULL
AFA_orginal_R2$gene_id <- as.character(AFA_orginal_R2$gene_id)

gencodev18 <- read.table(file = "/home/paul/mesa_models/gencode.v18.annotation.parsed.txt", header = T)
gencodev18$gene_id <- as.character(gencodev18$gene_id)

oldafa22 <- inner_join(AFA_orginal_R2, gencodev18, by = c("gene_id" = "gene_id"))
oldafa22 <- subset(oldafa22, oldafa22$chr == 22)

oldafa21 <- inner_join(AFA_orginal_R2, gencodev18, by = c("gene_id" = "gene_id"))
oldafa21 <- subset(oldafa21, oldafa21$chr == 21)

oldafa20 <- inner_join(AFA_orginal_R2, gencodev18, by = c("gene_id" = "gene_id"))
oldafa20 <- subset(oldafa20, oldafa20$chr == 20)

oldafa19 <- inner_join(AFA_orginal_R2, gencodev18, by = c("gene_id" = "gene_id"))
oldafa19 <- subset(oldafa19, oldafa19$chr == 19)

oldafa18 <- inner_join(AFA_orginal_R2, gencodev18, by = c("gene_id" = "gene_id"))
oldafa18 <- subset(oldafa18, oldafa18$chr == 18)

elnet_col <- colnames(afa_chr13)

#Filter the afa_chr13 and 14 by R2
afa_chr13 <- subset(afa_chr13, afa_chr13$cv_R2_avg >= 0.01)
afa_chr13$gene_name <- as.character(afa_chr13$gene_name)

afa_chr14 <- subset(afa_chr14, afa_chr14$cv_R2_avg >= 0.01)
afa_chr14$gene_name <- as.character(afa_chr14$gene_name)

#Figure out the overlap between prot_chr and afa_chr
afa_chr13_overlap <- semi_join(afa_chr13, prot_afa_chr13, by = c("gene_name" = "gene_name"))
afa_chr14_overlap <- semi_join(afa_chr14, prot_afa_chr14, by = c("gene_name" = "gene_name"))
prot_chr14_overlap <- semi_join(prot_afa_chr14, afa_chr14, by = c("gene_name" = "gene_name"))

cor.test(prot_afa_chr13$test_R2_avg, afa_chr13_overlap$test_R2_avg)
plot(prot_afa_chr13$test_R2_avg, afa_chr13_overlap$test_R2_avg, xlim = c(-0.5,1), ylim = c(-0.5,1), main = "AFA CHR13 models (PC_2,3 Vs PC_1,2,3)",
     xlab = "test_R2_avg chr13 PC_2,3", ylab = "test_R2_avg chr13 PC_1,2,3", frame=FALSE)
#abline(lm(afa_chr13_overlap$test_R2_avg ~ prot_afa_chr13$test_R2_avg), col = "blue")
cor.test(afa_chr14_overlap$test_R2_avg, prot_chr14_overlap$test_R2_avg)
plot(prot_chr14_overlap$test_R2_avg, afa_chr14_overlap$test_R2_avg, xlim = c(-0.5,1), ylim = c(-0.5,1), main = "AFA CHR14 models (PC_2,3 Vs PC_1,2,3)",
     xlab = "test_R2_avg chr14 PC_2,3", ylab = "test_R2_avg chr14 PC_1,2,3", frame=FALSE)

#histogram distribution of prot13 and 14, and afa 13 and 14
hist(prot_afa_chr13$test_R2_avg, breaks = "FD", xlim = c(-1, 1),freq = T, xlab = "R^2", main = "AFA Chr13 PC_2,3 test R2")
hist(prot_afa_chr14$test_R2_avg, breaks = "FD", xlim = c(-1, 1),freq = T, xlab = "R^2", main = "AFA Chr14 PC_2,3 test R2")

hist(afa_chr13$test_R2_avg, breaks = "FD", xlim = c(-1, 1),freq = T, xlab = "R^2", main = "AFA Chr13 PC_1,2,3 test R2")
hist(afa_chr14$test_R2_avg, breaks = "FD", xlim = c(-1, 1),freq = T, xlab = "R^2", main = "AFA Chr14 PC_1,2,3 test R2")

#Test the new AFA model window = 100, @ =0.05, seed = 1000, protein_coding only
afa_chr22 <- read.table(file = "/home/paul/mesa_models/train_results/new_AFA_model_nested_cv_chr22_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr21 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr21_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr20 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr20_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr19 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr19_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr18 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr18_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr14 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr14_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr15 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr15_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr16 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr16_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr17 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr17_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr13 <- read.table(file = "Z:/mesa_models/train_results/new_AFA_model_nested_cv_chr13_model_summaries_10_peer_3pcs.txt", header = TRUE)

afa_chr22 <- read.table(file = "/home/paul/mesa_models/train_results/new_AFA_model_nested_cv_chr22_model_summaries_10_peer_3pcs.txt", header = TRUE)
afa_chr22_weights <- read.table(file = "/home/paul/mesa_models/train_results/new_AFA_model_nested_cv_chr22_weights_10_peer_3pcs.txt", header = TRUE)

#Splitting the gene expression of every chromosome and having ryan to do same on the genotype
library(dplyr) 
# or
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")

gen18 <- read.table(file="Z:/mesa_models/gencode.v18.annotation.parsed.txt", header = T)

gex <- read.table(file="Z:/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
gex$PROBE_ID <- as.character(gex$PROBE_ID)

chr22_gencode <- subset(gen18, gen18$chr == 13)
chr22_gencode <- subset(chr22_gencode, chr22_gencode$gene_type == "protein_coding")
chr22_gencode$gene_id <- as.character(chr22_gencode$gene_id)
chr22_ex <- inner_join(chr22_gencode, gex, by = c("gene_id" = "PROBE_ID"))

#afa_chr12 <- read.table(file = "/home/paul/mesa_models/train_results/new_AFA_model_nested_cv_chr12_model_summaries_10_peer_3pcs.txt", header = TRUE)
#mydata[order(mydata$B),] how to order a dataframe by a column

chr22_genes <- chr22_ex
chr22_genes <- chr22_genes[order(chr22_genes$start),]

chr22_chunk1 <- chr22_genes[1:49,]
chr22_chunk1 <- chr22_chunk1[,-c(1,3:6)]
#chunk_col <- colnames(chr14_chunk1)
write.table(chr22_chunk1, file = "/home/paul/mesa_models/split_mesa/AFA_chr22_gex_chunk1.txt", row.names = FALSE)

chr22_chunk2 <- chr22_genes[50:98,]
chr22_chunk2 <- chr22_chunk2[,-c(1,3:6)]
write.table(chr22_chunk2, file = "/home/paul/mesa_models/split_mesa/AFA_chr22_gex_chunk2.txt", row.names = FALSE)

chr22_chunk3 <- chr22_genes[99:147,]
chr22_chunk3 <- chr22_chunk3[,-c(1,3:6)]
write.table(chr22_chunk3, file = "/home/paul/mesa_models/split_mesa/AFA_chr22_gex_chunk3.txt", row.names = FALSE)

chr22_chunk4 <- chr22_genes[148:196,]
chr22_chunk4 <- chr22_chunk4[,-c(1,3:6)]
write.table(chr22_chunk4, file = "/home/paul/mesa_models/split_mesa/AFA_chr22_gex_chunk4.txt", row.names = FALSE)

chr22_chunk5 <- chr22_genes[197:243,]
chr22_chunk5 <- chr22_chunk5[,-c(1,3:6)]
write.table(chr22_chunk5, file = "/home/paul/mesa_models/split_mesa/AFA_chr22_gex_chunk5.txt", row.names = FALSE)



#AFA Gene Expression Splitting loop
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")

gen18 <- read.table(file="/home/paul/mesa_models/gencode.v18.annotation.parsed.txt", header = T)

gex <- read.table(file="/home/paul/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
gex$PROBE_ID <- as.character(gex$PROBE_ID)
pop <- "AFA"
for (i in 1:22){
  chr_gencode <- subset(gen18, gen18$chr == i)
  chr_gencode <- subset(chr_gencode, chr_gencode$gene_type == "protein_coding")
  chr_gencode$gene_id <- as.character(chr_gencode$gene_id)
  chr_ex <- inner_join(chr_gencode, gex, by = c("gene_id" = "PROBE_ID"))
  chr_genes <- chr_ex
  chr_genes <- chr_genes[order(chr_genes$start),]
  n_genes <- nrow(chr_genes)
  n <- floor(n_genes/5)
  ch <- n*5 #recalculate the nrow
  d <- n_genes - ch #take the remainder
  l <- n + (d) #add the remainder to be the last number
  va <- c(n,n,n,n,l)# just vectors of the division numbers
  vi <- c(1, n+1, n+1+n, n+1+n+n, n+1+n+n+n) #vector of starting values
  vf <- c(n, n+n, n+n+n, n+n+n+n, n+n+n+n+l) #vectors of ending values
  i <- as.character(i)
  for (j in 1:5){
    chr_chunk <- chr_genes[vi[j]:vf[j],]
    chr_chunk <- chr_chunk[,-c(1,3:6)]
    j <- as.character(j)
    write.table(chr_chunk, file = "/home/paul/mesa_models/split_mesa/" %&% pop %&% "_chr" %&% i %&% "_gex_chunk" %&% j %&% ".txt", row.names = FALSE, quote = FALSE, sep = "\t")
  }
}





#CAU Gene Expression Splitting
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")

gen18 <- read.table(file="/home/paul/mesa_models/gencode.v18.annotation.parsed.txt", header = T)

gex <- read.table(file="/home/paul/mesa_models/meqtl_sorted_CAU_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
gex$PROBE_ID <- as.character(gex$PROBE_ID)
pop <- "CAU"
for (i in 1:22){
  chr_gencode <- subset(gen18, gen18$chr == i)
  chr_gencode <- subset(chr_gencode, chr_gencode$gene_type == "protein_coding")
  chr_gencode$gene_id <- as.character(chr_gencode$gene_id)
  chr_ex <- inner_join(chr_gencode, gex, by = c("gene_id" = "PROBE_ID"))
  chr_genes <- chr_ex
  chr_genes <- chr_genes[order(chr_genes$start),]
  n_genes <- nrow(chr_genes)
  n <- floor(n_genes/5)
  ch <- n*5 #recalculate the nrow
  d <- n_genes - ch #take the remainder
  l <- n + (d) #add the remainder to be the last number
  va <- c(n,n,n,n,l)# just vectors of the division numbers
  vi <- c(1, n+1, n+1+n, n+1+n+n, n+1+n+n+n) #vector of starting values
  vf <- c(n, n+n, n+n+n, n+n+n+n, n+n+n+n+l) #vectors of ending values
  i <- as.character(i)
  for (j in 1:5){
    chr_chunk <- chr_genes[vi[j]:vf[j],]
    chr_chunk <- chr_chunk[,-c(1,3:6)]
    j <- as.character(j)
    write.table(chr_chunk, file = "/home/paul/mesa_models/split_mesa/" %&% pop %&% "_chr" %&% i %&% "_gex_chunk" %&% j %&% ".txt", row.names = FALSE, quote = FALSE, sep = "\t")
  }
}







#HIS Gene Expression Splitting
library(tidyverse)
"%&%" <- function(a,b) paste(a,b, sep = "")

gen18 <- read.table(file="/home/paul/mesa_models/gencode.v18.annotation.parsed.txt", header = T)

gex <- read.table(file="/home/paul/mesa_models/meqtl_sorted_HIS_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
gex$PROBE_ID <- as.character(gex$PROBE_ID)
pop <- "HIS"
for (i in 1:22){
  chr_gencode <- subset(gen18, gen18$chr == i)
  chr_gencode <- subset(chr_gencode, chr_gencode$gene_type == "protein_coding")
  chr_gencode$gene_id <- as.character(chr_gencode$gene_id)
  chr_ex <- inner_join(chr_gencode, gex, by = c("gene_id" = "PROBE_ID"))
  chr_genes <- chr_ex
  chr_genes <- chr_genes[order(chr_genes$start),]
  n_genes <- nrow(chr_genes)
  n <- floor(n_genes/5)
  ch <- n*5 #recalculate the nrow
  d <- n_genes - ch #take the remainder
  l <- n + (d) #add the remainder to be the last number
  va <- c(n,n,n,n,l)# just vectors of the division numbers
  vi <- c(1, n+1, n+1+n, n+1+n+n, n+1+n+n+n) #vector of starting values
  vf <- c(n, n+n, n+n+n, n+n+n+n, n+n+n+n+l) #vectors of ending values
  i <- as.character(i)
  for (j in 1:5){
    chr_chunk <- chr_genes[vi[j]:vf[j],]
    chr_chunk <- chr_chunk[,-c(1,3:6)]
    j <- as.character(j)
    write.table(chr_chunk, file = "/home/paul/mesa_models/split_mesa/" %&% pop %&% "_chr" %&% i %&% "_gex_chunk" %&% j %&% ".txt", row.names = FALSE, quote = FALSE, sep = "\t")
  }
}



n_genes <- nrow(chr22_ex)
n_genes <- 1048
n <- floor(n_genes/5)
ch <- n*5
d <- n_genes - ch
l <- n + (d)
va <- c(n,n,n,n,l)

#Test the spliited gex files and see they are correct
spl_test <- read.table(file = "Z:/mesa_models/split_mesa/AFA_chr2_gex_chunk1.txt", header = T)



#Test the chunk results
chr1chunk1 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr1_model_summaries_chunk1.txt", header = TRUE)
chr1chunk2 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr1_model_summaries_chunk2.txt", header = TRUE)
chr1chunk3 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr1_model_summaries_chunk3.txt", header = TRUE)
chr1chunk4 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr1_model_summaries_chunk4.txt", header = TRUE)
chr1chunk5 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr1_model_summaries_chunk5.txt", header = TRUE)
full_AFA_chr13 <- read.table(file = "Z:/mesa_models/split_mesa/results/full_AFA_chr13_model_summaries.txt", header = TRUE)
full_AFA_chr14 <- read.table(file = "Z:/mesa_models/split_mesa/results/full_AFA_chr14_model_summaries.txt", header = TRUE)
full_AFA_chr16 <- read.table(file = "Z:/mesa_models/split_mesa/results/full_AFA_chr16_model_summaries.txt", header = TRUE)


##See if these genes are expressed and with positive R2 in MESA: FFAR1, FFAR2, FFAR3, FFAR4
AFA_orginal_R2 <- read.table(file = "/home/paul/mesa_models/AFA_model_summary_stats.txt", header = F)
AFA_chr10_elnet <- read.table(file = "/home/paul/mesa_models/uncompleted_old_train_results/new_AFA_model_nested_cv_chr10_model_summaries_10_peer_3pcs.txt", header = T)


chr1chunk2 <- read.table(file = "/home/paul//mesa_models/split_mesa/results/split_AFA_chr1_model_summaries_chunk2.txt", header = TRUE)

chr11chunk1 <- read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr11_model_summaries_chunk1.txt", header = TRUE)
chr11chunk1we <- read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr11_weights_chunk1.txt", header = TRUE)
chr11chunk1cov <- read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr11_covariances_chunk1.txt", header = TRUE)
chr11chunk2cov <- read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr11_covariances_chunk2.txt", header = TRUE)
chr11chunk3cov <- read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr11_covariances_chunk3.txt", header = TRUE)
chr11chunk4cov <- read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr11_covariances_chunk4.txt", header = TRUE)


chr10chunk2 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr10_model_summaries_chunk2.txt", header = TRUE)
chr10chunk3 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr10_model_summaries_chunk3.txt", header = TRUE)
chr10chunk4 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr10_model_summaries_chunk4.txt", header = TRUE)
chr10chunk5 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr10_model_summaries_chunk5.txt", header = TRUE)

chr14chunk1 <- read.table(file = "Z:/mesa_models/split_mesa/results/split_AFA_chr14_model_summaries_chunk1.txt", header = TRUE)


#combine the doMC mesa chunk results
"%&%" <- function(a,b) paste(a,b, sep='')
for (chrom in c(1:12, 14, 16, 17, 19)){
  #summaries <- NULL
  #weights <- NULL
  covariances <- NULL
  for (chunk in 1:5){
    #summaries <- rbind(summaries, read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr" %&% as.character(chrom) %&% "_model_summaries_chunk" %&% as.character(chunk) %&% ".txt", header = T, stringsAsFactors = F))
    #weights <- rbind(weights, read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr" %&% as.character(chrom) %&% "_weights_chunk" %&% as.character(chunk) %&% ".txt", header = T, stringsAsFactors = F))
    covariances <- rbind(covariances, read.table(file = "/home/paul/mesa_models/split_mesa/results/split_AFA_chr" %&% as.character(chrom) %&% "_covariances_chunk" %&% as.character(chunk) %&% ".txt", header = T, stringsAsFactors = F))
  }
  #write.table(summaries, file = "/home/paul/mesa_models/split_mesa/results/full_AFA_chr" %&% as.character(chrom) %&% "_model_summaries.txt", row.names = FALSE)
  #write.table(weights, file = "/home/paul/mesa_models/split_mesa/results/full_AFA_chr" %&% as.character(chrom) %&% "_weights.txt", row.names = FALSE)
  write.table(covariances, file = "/home/paul/mesa_models/split_mesa/results/full_AFA_chr" %&% as.character(chrom) %&% "_covariances.txt", row.names = FALSE)
}


#join all chromosomes to create one file for model results
model_summaries <- NULL
for (chrom in 1:22) {
  model_summaries <- rbind(model_summaries, read.table(file = "/home/paul/mesa_models/split_mesa/results/full_AFA_chr" %&% as.character(chrom) %&% "_model_summaries.txt", header = T, stringsAsFactors = F))
}

#Just write out the model_summaries
write.table(model_summaries, file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", row.names = FALSE)

#filter model summaries by cv R2
model_summaries <- subset(model_summaries, model_summaries$cv_R2_avg > 0.5)

mesa_all <- read.table(file = "Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)

hist(mesa_all$cv_R2_avg, breaks = "FD", xlim = c(-1, 1), freq = T, xlab = "CV R2", col = "green", main = "MESA AFA All Chromosomes Models Performances")

mesa_all <- subset(mesa_all,cv_R2_avg > 0.1)

#HIS model
mesa_his <- read.table(file = "/home/ryan/enet_scripts/HIS_MESA_models/MESA.HIS.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header = TRUE)
mesa_his <- subset(mesa_his, cv_R2_avg > 0.1)

mesa_cau <- read.table(file = "/home/ryan/enet_scripts/CAU_MESA_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header = TRUE)
mesa_cau <- subset(mesa_cau, cv_R2_avg > 0.1)

mesa2mets <- read.table(file = "/home/paul/pearson_vs_R2.txt", header = TRUE)

mesa2mets <- subset(mesa2mets, pearson > 0.1)
#extract the genenames from chr22
mesachr17 <- read.table(file = "/home/paul/mesa_models/split_mesa/results/full_AFA_chr17_model_summaries.txt", header = TRUE)
mesachr16 <- read.table(file = "/home/paul/mesa_models/split_mesa/results/full_AFA_chr16_model_summaries.txt", header = TRUE)

genelist <- as.character(mesachr22$gene_name)
#write out the gene names to a text file 
fileConn<-file("/home/paul/mesa_models/svr/genelist.txt")
writeLines(genelist, fileConn)
close(fileConn)

gg <- read.csv(file = "/home/paul/mesa_models/svr/svr22/svr_cis_gt_chr2_OTOF.csv") #Just to view the cis_gt

#extract the genenames from each of chr1 - 21
"%&%" <- function(a,b) paste(a,b, sep = "")
for (chrom in 1:21) {
  model <- read.table(file = "/home/paul/mesa_models/split_mesa/results/full_AFA_chr" %&% as.character(chrom) %&% "_model_summaries.txt", header = T)
  genelist <- as.character(model$gene_name)
  fileConn<-file("/home/paul/mesa_models/svr/genelist_chr" %&% as.character(chrom) %&% ".txt")
  writeLines(genelist, fileConn)
  close(fileConn)
}
#extract the genenames and r2 of the glmnet model 
relnet <- mesachr18[,c(2,10)]
#read in the python elnet
py_elnet18 <- read.table(file = "/home/paul/mesa_models/svr/chr18_CV_R2.txt", header = T)

#sort both the glmnet and pyelnet dataframes by gene names
#df <- df[order(df$x),]
relnet <- relnet[order(relnet$gene_name),]
py_elnet18 <- py_elnet18[order(py_elnet18$Gene_Name),]

#plot the relnet and pyelnet R2 against each other
plot(relnet$cv_R2_avg, py_elnet22$CV_R2_avg, xlab = "R glmnet model", ylab = "python sklearn elnet", 
     main = "python elnet 5 fold CV Vs R glmnet Nested CV for chr22")
abline(0,1, col="blue") #+
abline(-0.09009, 1.15280, col="red")

lm(r_py$CVP ~ r_py$CVR)
cor.test(relnet$cv_R2_avg, py_elnet22$CV_R2_avg)

#or

r_py <- cbind(relnet, py_elnet18)
colnames(r_py) <- c("gene_nameR", "CVR", "gene_nameP", "CVP")
library("ggpubr")
ggscatter(r_py, x = "CVR", y = "CVP", 
          add = "reg.line", add.params = list(color="blue"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "R glmnet model", ylab = "python sklearn elnet", title = "python elnet 5 fold CV Vs R glmnet Nested CV for chr18",
          xlim = c(-0.25, 0.75), ylim = c(-0.25, 0.75)) + geom_abline(intercept = 0, slope = 1, color="red")

p <- ggscatter(r_py, x = "CVR", y = "CVP", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "R glmnet model", ylab = "python sklearn elnet", title = "python elnet 5 fold CV Vs R glmnet Nested CV for chr21")

gpar(p, xlim = c(-0.25, 0.75)) #use ggpar to edit the parameters such as xlim and ylim
p + geom_abline(intercept = 0, slope = 1)


#expression files are from lauren
#/media/MyBook/lauren/files_for_revisions_plosgen/en_v7/prepare_data/expression


CAU <- read.table(file = "Z:/METS_model/METS_2_CAU_spearman2.txt", header = T)
ALL <- read.table(file = "Z:/METS_model/METS_2_ALL_spearman2.txt", header = T)
HIS <- read.table(file = "Z:/METS_model/METS_2_HIS_spearman2.txt", header = T)
AFA <- read.table(file = "Z:/METS_model/METS_2_AFA_spearman2.txt", header = T)
AFHI <- read.table(file = "Z:/METS_model/METS_2_AFHI_spearman2.txt", header = T)

library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

# colnames(ALL)<-c("gene.ALL","spearman.ALL")
# colnames(AFHI)<-c("gene.AFHI","spearman.AFHI")
# colnames(AFA)<-c("gene.AFA","spearman.AFA")
# colnames(CAU)<-c("gene.CAU","spearman.CAU")
# colnames(HIS)<-c("gene.HIS","spearman.HIS")

ALL$pop<-"ALL"
ALL$median<-round(median(ALL$spearman),3)
AFHI$pop<-"AFHI"
AFHI$median<-round(median(AFHI$spearman),3)
AFA$pop<-"AFA"
AFA$median<-round(median(AFA$spearman),3)
CAU$pop<-"CAU"
CAU$median<-round(median(CAU$spearman),3)
HIS$pop<-"HIS"
HIS$median<-round(median(HIS$spearman),3)
total<-rbind.data.frame(ALL,AFHI,AFA,CAU,HIS)
total[is.na(total)]<-0
pV<-ggplot(data=total,aes(y=spearman,x=as.factor(median)))
pV + geom_violin(draw_quantiles = T, aes(fill=pop)) + 
  geom_boxplot(width = 0.4) +
  geom_hline(yintercept=0) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("MESA spearman correlation using METS models chr 1, 9 - 22")


caugex_ryan <- read.table(file = "/home/paul/METS_model/CAU_PF10.txt", header = T)
caugex_paul <- read.table(file="/home/paul/mesa_models/meqtl_sorted_CAU_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)

pcoryp <- data.frame(gene = NULL, pearson = NULL) #create empty dataframe to store results
for (i in 1:length(caugex_ryan$PROBE_ID)){
  pcoryp[i,1] <- caugex_ryan$PROBE_ID[i]
  pcoryp[i,2] <- cor(as.numeric(caugex_paul[i,2:579]), as.numeric(caugex_ryan[i,2:579]), method = "pearson")
}
colnames(pcoryp) <- c("gene", "pearson")

plot(as.numeric(caugex_paul[1,2:579]), as.numeric(caugex_ryan[1,2:579]))

afagex_ryan <- read.table(file = "/home/paul/METS_model/AFA_PF10.txt", header = T)
afagex_paul <- read.table(file="/home/paul/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
cor(as.numeric(afagex_paul[1,2:234]), as.numeric(afagex_ryan[1,2:234]))

#check if the columns are sorted
#test[ , order(names(test))]

#df1 %>% dplyr::select(vector_of_column_names_you_want_to_match)
afagex_paul <- afagex_paul %>% dplyr::select(names(afagex_ryan)) #make the colnames to be in matching order
sum(colnames(afagex_paul) == colnames(afagex_ryan)) #check if they are in the same matching order

caugex_paul <- caugex_paul[, order(names(caugex_paul))] #another way of sorting columns
caugex_ryan <- caugex_ryan[, order(names(caugex_ryan))] #another way of sorting columns
sum(colnames(caugex_paul) == colnames(caugex_ryan)) #confirm again


mesa_all <- read.table(file = "/home/paul/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header = TRUE)
knn <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/knn_cv_full_AFA_chr.txt", header = TRUE)
rf <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/rf_cv_full_AFA_chr.txt", header = TRUE)
svr_rbf <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/svr_rbf_cv_full_AFA_chr.txt", header = TRUE)
svr_linear <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/svr_linear_cv_full_AFA_chr.txt", header = TRUE)


split1 <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/grid_split/best_grid_split_knn_cv_chr10_chunk1.txt", header = T)
