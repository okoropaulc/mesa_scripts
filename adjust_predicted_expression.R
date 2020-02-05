#Adjust the predicted expression with 10 PCs
library(data.table)

#Read in the PCA eigenvector
pca10 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/unmerged_pca.eigenvec", header=T)
pca10$FID <- NULL

#make the ID to be numbers only
for (i in 1:nrow(pca10)){
  pca10$IID[i] <- strsplit(pca10$IID[i], "_")[[1]][2]
}
pca10$IID <- as.character(pca10$IID)

#keep only samples that have rankplt5
#read in pheno so as to know sample id to retain
thrombomodulin <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)
library(tidyverse)

#drop NA
#keep ID and rankplt5
thrombomodulin <- thrombomodulin[,c(1,5)]
thrombomodulin <- drop_na(thrombomodulin) #drop NA
thrombomodulin$sidno <- as.character(thrombomodulin$sidno)
colnames(thrombomodulin)[2] <- "phenotype" #rank_plt5

library(dplyr)
pcsid <- inner_join(thrombomodulin, pca10, by = c("sidno" = "IID"))
thrombomodulin <- pcsid[,c(1,2)] #store back only samples with pcs
cov_df <- pcsid[,(3:12)] #store only the 10 PCs. They are in the same sample order with the thrombomodulin because of the inner_join


#read in pred expr and retain only samples that are in the phenotype
predrf <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_expr.txt",header=T)
predsvr <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_svr_pred_expr.txt",header=T)
predknn <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_knn_pred_expr.txt",header=T)
pxcan <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/en_cau_rankplt5_predicted_expression.txt", header=T)
pxcan$FID <- NULL #remove the first column

#take the genes
pxcangenes <- colnames(pxcan)[2:length(pxcan)]
rfgenes <- colnames(predrf)[2:length(predrf)]
svrgenes <- colnames(predsvr)[2:length(predsvr)]
knngenes <- colnames(predknn)[2:length(predknn)]

predrf$sampleidrf <- as.character(predrf$sampleidrf)
predsvr$sampleidsvr <- as.character(predsvr$sampleidsvr)
predknn$sampleidknn <- as.character(predknn$sampleidknn)
pxcan$IID <- as.character(pxcan$IID)

library(dplyr)

#reduce the pred ex to have only sampe in the thrombomodulin and remove the pheno column
pxcan2 <- inner_join(thrombomodulin, pxcan, by = c("sidno" = "IID"))
pxcan2$phenotype <- NULL
rf2 <- inner_join(thrombomodulin, predrf, by = c("sidno" = "sampleidrf"))
rf2$phenotype <- NULL
svr2 <- inner_join(thrombomodulin, predsvr, by = c("sidno" = "sampleidsvr"))
svr2$phenotype <- NULL
knn2 <- inner_join(thrombomodulin, predknn, by = c("sidno" = "sampleidknn"))
knn2$phenotype <- NULL

#adjust predicted expressions with 10 pcs

adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

pxcan_adj <- NULL
for (gene in pxcangenes){
  gene_expr <- pxcan2[[gene]]
  pxcan_adj <- cbind(pxcan_adj,adjust_for_covariates(gene_expr, cov_df))
}
pxcan_adj <- as.data.frame(pxcan_adj)
names(pxcan_adj) <- pxcangenes

rf_adj <- NULL
for (gene in rfgenes){
  gene_expr <- rf2[[gene]]
  rf_adj <- cbind(rf_adj,adjust_for_covariates(gene_expr, cov_df))
}
rf_adj <- as.data.frame(rf_adj)
names(rf_adj) <- rfgenes

svr_adj <- NULL
for (gene in svrgenes){
  gene_expr <- svr2[[gene]]
  svr_adj <- cbind(svr_adj, adjust_for_covariates(gene_expr, cov_df))
}
svr_adj <- as.data.frame(svr_adj)
names(svr_adj) <- svrgenes

knn_adj <- NULL
for (gene in knngenes){
  gene_expr <- knn2[[gene]]
  knn_adj <- cbind(knn_adj, adjust_for_covariates(gene_expr, cov_df))
}
knn_adj <- as.data.frame(knn_adj)
names(knn_adj) <- knngenes

#merge pheno and PC adjusted pred expr on their common sample. write out the adjusted pred expr, do pheno assoc down stream
knn_adj <- cbind(thrombomodulin,knn_adj)
#some gene columns are empty, so find them and remove them
for (gene in knngenes){
  if (is.null(knn_adj[[gene]])){
    knn_adj$gene <- NULL
  }
}

wknn <- knn_adj[,c(1,3:length(knn_adj))] # keep id and genes #w = for writing out. the 9625 is the df ncol
fwrite(wknn, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/knn_PC10_adjusted_pred_expr.txt", quote=F, row.names=F, sep="\t")
aknn <- knn_adj[,c(2:length(knn_adj))] #keep pheno and genes #a = for assoc

svr_adj <- cbind(thrombomodulin,svr_adj)
wsvr <- svr_adj[,c(1,3:length(svr_adj))] # keep id and genes #w = for writing out
fwrite(wsvr, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/svr_PC10_adjusted_pred_expr.txt", quote=F, row.names=F, sep="\t")
asvr <- svr_adj[,c(2:length(svr_adj))] #keep pheno and genes #a = for assoc

rf_adj <- cbind(thrombomodulin,rf_adj)
wrf <- rf_adj[,c(1,3:length(rf_adj))] # keep id and genes #w = for writing out
fwrite(wrf, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_adjusted_pred_expr.txt", quote=F, row.names=F, sep="\t")
arf <- rf_adj[,c(2:length(rf_adj))] #keep pheno and genes #a = for assoc

pxcan_adj <- cbind(thrombomodulin,pxcan_adj)
wpx <- pxcan_adj[,c(1,3:length(pxcan_adj))] # keep id and genes #w = for writing out
fwrite(wpx, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_adjusted_pred_expr.txt", quote=F, row.names=F, sep="\t")
apx <- pxcan_adj[,c(2:length(pxcan_adj))] #keep pheno and genes #a = for assoc

# trypx <- apx
# for (gene in pxcangenes){
#   if (is.null(trypx[[gene]])){
#     trypx$gene <- NULL
#   }
# }

# functionize

#ascociation function

association <- function(merged, genes, test_type = "linear") {
  assoc_df <- NULL # Init association dataframe
  
  # Perform test between each pred_gene_exp column and phenotype
  for (gene in genes) {
    pred_gene_exp <- merged[[gene]]
    if (test_type == "logistic") { 
      model <- glm(phenotype ~ pred_gene_exp, data = merged, family = binomial)
    } else if (test_type == "linear") {
      model <- lm(phenotype ~ pred_gene_exp, data = merged)
    } else if (test_type == "survival") {
      # TODO: survival analysis
      model <- NULL
    }
    results <- coef(summary(model))[c(2,6,8,4)]
    line <- c(gene,results)
    assoc_df <- rbind(assoc_df,line)
  }
  
  # Specify column names of assoc_df
  if (test_type == "logistic") {
    colnames(assoc_df) <- c("gene", "beta", "z", "p", "se(beta)")
  } else if (test_type == "linear") {
    colnames(assoc_df) <- c("gene", "beta", "t", "p", "se(beta)")
  } else if (test_type == "survival") {
    # TODO
  }
  return(as.data.frame(assoc_df))
}

pxcanassoc <- association(apx, pxcangenes)
rfassoc <- association(arf, rfgenes)
svrassoc <- association(asvr, svrgenes)
knnassoc <- association(aknn, knngenes)

#write out the result
fwrite(rfassoc, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_rankplt5_association.txt",row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/svr_PC10_rankplt5_association.txt",row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/knn_PC10_rankplt5_association.txt",row.names=F, quote=F, sep="\t")
fwrite(pxcanassoc, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_rankplt5_association.txt", quote=F, row.names=F, sep ="\t")


#make manhattan plot from the associations. also filter out genes with cv_r2 < 0.01 for the ML

library(data.table)
library(dplyr)
#note that there is a bug in fread which makes it to skip some row for huge files

rfassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_rankplt5_association.txt", header=T)
svrassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/svr_PC10_rankplt5_association.txt", header=T)
knnassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/knn_PC10_rankplt5_association.txt", header=T)


#make qqplot
#install.packages("qqman")
library(qqman)

qq(knnassoc$p, main="knn PC10 unfiltered")
qq(svrassoc$p, main="svr PC10 unfiltered")
qq(rfassoc$p, main="rf PC10 unfiltered")

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
fwrite(rfassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/svr_PC10_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/knn_PC10_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")


#ggplot after filtering
qq(rfassoc_cvr0.01$p, main="rf PC10 filtered by cv r2")
qq(svrassoc_cvr0.01$p, main="svr PC10 filtered by cv r2")
qq(knnassoc_cvr0.01$p, main="knn PC10 filtered by cv r2")
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
rfassocgen <- rfassocgen[,c(1,3,4,9)]#drop gene_id and use gene_name to help annotation in the plot
rfassocgen <- inner_join(gencode1, rfassoc, by = c("gene_id" = "gene")) #unfiltered
rfassocgen <- rfassocgen[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svrassocgen <- inner_join(gencode1, svrassoc_cvr0.01, by = c("gene_id" = "Gene_ID"))
svrassocgen <- svrassocgen[,c(1,3,4,9)] #drop gene_id
svrassocgen <- inner_join(gencode1, svrassoc, by = c("gene_id" = "gene"))
svrassocgen <- svrassocgen[,c(1,3,4,7)] #drop gene_id

knnassocgen <- inner_join(gencode1, knnassoc_cvr0.01, by = c("gene_id" = "Gene_ID"))
knnassocgen <- knnassocgen[,c(1,3,4,9)] #drop gene_id
knnassocgen <- inner_join(gencode1, knnassoc, by = c("gene_id" = "gene")) #manhattan for unfiltered
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
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="rf PC10 adjusted No cv_R2 filtering", annotatePval=0.0005,
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
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="svr PC10 No cv_R2 filtering", annotatePval=0.0005,
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
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="knn PC10 No cv_R2 filtering", 
          annotatePval=0.0005, col=c("red", "blue"))




#Predixcan
enassoc <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_rankplt5_association.txt", header=T)
enassoc$gene <- as.character(enassoc$gene)
for (i in 1:nrow(enassoc)){
  enassoc$gene[i] <- gsub('\\.[0-9]+','',enassoc$gene[i])
} #just to remove the decimal places in the gene_id


#make qqplot
#install.packages("qqman")
library(qqman)

qq(enassoc$p, main="en PC10 unfiltered")

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
fwrite(enassoc_cvr0.01, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_assoc_filtered_by_r2_0.01.txt", row.names=F, quote=F, sep="\t")

#ggplot after filtering
qq(enassoc_cvr0.01$p, main="en PC10 filtered by cv r2")


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

enassocgen <- inner_join(gencode1, enassoc_cvr0.01, by = c("gene_id" = "gene_id"))
enassocgen <- enassocgen[,c(1,3,4,9)] #drop gene_id and use gene_name to help annotation in the plot
enassocgen <- inner_join(gencode1, enassoc, by = c("gene_id" = "gene")) #do unfiltered
enassocgen <- enassocgen[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

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
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="predixcan PC10 No cv_R2 filtering corrected", annotatePval=0.0005,
          col=c("red", "blue"))





###################################################
########################################################################################
####################################################################################################################

#check the direction of effect for the hit genes in rf(TRPM4=ENSG00000130529.11) and en and svr(MCM3AP=ENSG00000160294.6)
#knn (TMEM50B=ENSG00000142188.12)
library(data.table)
#check gene name
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
# rf(ZFP36L1=ENSG00000185650.8, sox15=ENSG00000129194.3). and en (DESI1=ENSG00000100418.7)
#read in the prediction from both algs

pxcan <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_adjusted_pred_expr.txt", header=T)
predrf <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_adjusted_pred_expr.txt",header=T)
predsvr <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/svr_PC10_adjusted_pred_expr.txt",header=T)
predknn <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/knn_PC10_adjusted_pred_expr.txt",header=T)

#extract the predictd expression for the genes
en_mcm3ap <- pxcan[["ENSG00000160294.6"]]
en_desi1 <- pxcan[["ENSG00000100418.7"]]
rf_trpm4 <-predrf[["ENSG00000130529.11"]]
rf_zfp36l1 <- predrf[["ENSG00000185650.8"]]
rf_sox15 <- predrf[["ENSG00000129194.3"]]
svr_mcm3ap <- predsvr[["ENSG00000160294"]] #mcm3ap gene id does not have . in svr
knn_tmem50b <- predknn[["ENSG00000142188"]] #tmem50b gene id does not have . in knn
pheno <- thrombomodulin$phenotype #the phenotype is same for all merged df

#make a df and store pheno and genes of interest expression values. it will be easy to plot them that way
pheno_g <- data.frame(pheno=pheno,en_mcm3ap=en_mcm3ap,en_desi1=en_desi1,rf_trpm4=rf_trpm4,rf_zfp36l1=rf_zfp36l1,rf_sox15=rf_sox15,
                      svr_mcm3ap=svr_mcm3ap, knn_tmem50b=knn_tmem50b)

#plot them
#MCM3AP en
library(ggplot2)
ggplot(pheno_g, aes(x=en_mcm3ap, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") +
  theme_classic(20) + ggtitle("Elastic Net Gene MCM3AP PC10 Adjusted") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) +
  geom_density_2d(lwd=1.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#MCM3AP svr
library(ggplot2)
ggplot(pheno_g, aes(x=svr_mcm3ap, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") +
  theme_classic(20) + ggtitle("SVR Gene MCM3AP PC10 Adjusted") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) +
  geom_density_2d(lwd=1.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#DESI1
library(ggplot2)
ggplot(pheno_g, aes(x=en_desi1, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") +
  theme_classic(20) + ggtitle("Elastic Net Gene DESI1 PC10 Adjusted") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#TRPM4
library(ggplot2)
ggplot(pheno_g, aes(x=rf_trpm4, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") +
  theme_classic(20) + ggtitle("Random Forest Gene TRPM4 PC10 Adjusted") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) +
  geom_density_2d(lwd=1.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#ZFP36L1
library(ggplot2)
ggplot(pheno_g, aes(x=rf_zfp36l1, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") +
  theme_classic(20) + ggtitle("Random Forest Gene ZFP36L1 PC10 Adjusted") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) +
  geom_density_2d(lwd=1.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#SOX15
library(ggplot2)
ggplot(pheno_g, aes(x=rf_sox15, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") +
  theme_classic(20) + ggtitle("Random Forest Gene SOX15 PC10 Adjusted") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

#TMEM50B
library(ggplot2)
ggplot(pheno_g, aes(x=knn_tmem50b, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("platelet count (rank normalized)") +
  theme_classic(20) + ggtitle("KNN Gene TMEM50B PC10 Adjusted") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")


# plot(rf_trpm4,pheno)
# plot(rf_zfp36l1, pheno)
# plot(rf_sox15, pheno)
# plot(en_mcm3ap, pheno)
# plot(en_desi1, pheno)

#compare the t statistic of rf and en on their overlap genes
# read in the filtered assoc files and use gene_name to overlap them
rfassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_assoc_filtered_by_r2_0.01.txt", header=T)
rfassoc_cvr0.01$Gene_Name <- as.character(rfassoc_cvr0.01$Gene_Name)
enassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_assoc_filtered_by_r2_0.01.txt", header=T)
enassoc_cvr0.01$gene_name <- as.character(enassoc_cvr0.01$gene_name)

rfen <- inner_join(rfassoc_cvr0.01, enassoc_cvr0.01, by =c("Gene_Name"="gene_name"))
rfen <- rfen[,c(5,11)] #take out their t statistic
names(rfen) <- c("rf","en")
# #first remove decimals from the gene id
# for (i in 1:length(rf_afa$Gene_ID)){
#   rf_afa$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_afa$Gene_ID[i])
# } #just to remove the decimal places in the gene_id

library(ggplot2)
ggplot(rfen, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("Random Forest") + xlim(-6,6) + ylim(-6,6) +
  theme_classic(20) + ggtitle("t-statistic comparison PC10 Adjusted") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)

#use ggpubr
library("ggpubr")
ggscatter(rfen, x = "en", y = "rf", 
          add = "reg.line", add.params = list(color="red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elastic Net", ylab = "Random Forest", 
          title = "t-statistic comparison PC10 Adjusted",
          xlim = c(-6, 6), ylim = c(-6, 6)) + geom_abline(intercept = 0, slope = 1, color="blue")



#read in the predicted expressions

library(data.table)

wknn <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/knn_PC10_adjusted_pred_expr.txt", header=T)
wsvr <- fread(wsvr, file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/svr_PC10_adjusted_pred_expr.txt", header=T)
wrf <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_adjusted_pred_expr.txt", header=T)
wpx <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_adjusted_pred_expr.txt", header=T)


enassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/en_PC10_assoc_filtered_by_r2_0.01.txt", header=T)
rfassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/rf_PC10_assoc_filtered_by_r2_0.01.txt", header=T)
svrassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/svr_PC10_assoc_filtered_by_r2_0.01.txt",header=T)
knnassoc_cvr0.01 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/adjusted_pred_expr/knn_PC10_assoc_filtered_by_r2_0.01.txt",header=T)

