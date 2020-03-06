#This script will use the pcs to adjust the predicted expression
#it will also run twas of the expression and phenotype

library(data.table)
library(dplyr)

#read in the pcs of the ld pruned mesa ALL
mesa_pc <- read.table(file="Z:/data/lauren_mesa/allpops_dosages_joined/ld_pruned_mesa_unmerged.eigenvec", header=T)

#read in the predicted expression already filtered by r2 of each algorithm
rf <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_rf.txt", header = T)
svr <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_svr.txt", header = T)
knn <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_knn.txt", header = T)
en <- fread(file = "Z:/data/twas_mesa/ALL_en_predicted_expression_r2_0.01.txt", header = T)

#keep copy of the files, just to keep a copy of these files because it takes long to read them in rstudio
rf_cp <-rf
svr_cp <- svr
knn_cp <- knn
en_cp <- en

#Now arrange the sample IID of the predicted expressions to be in same order with the mesa_pc sample IID. arrange pheno with this too
#use inner_join
IID_order <- data.frame(IID=mesa_pc[,2])
fwrite(IID_order, file="Z:/data/twas_mesa/samples_for_adj_pred_exp.txt", row.names=F, quote=F)

rf <- inner_join(IID_order, rf, by = c("IID"="IID"))
svr <- inner_join(IID_order, svr, by = c("IID"="IID"))
knn <- inner_join(IID_order, knn, by = c("IID"="IID"))
en <- inner_join(IID_order, en, by = c("IID"="IID"))

cov_df <- mesa_pc[,(3:12)]#take only the 10 pcs since the pc sample IID is in same order with the predicted expression

#take the genes names in the predicted expression dataframe, so I can use it to iteratively adjust the expression with pc
pxcangenes <- colnames(en)[2:length(en)]
rfgenes <- colnames(rf)[2:length(rf)]
svrgenes <- colnames(svr)[2:length(svr)]
knngenes <- colnames(knn)[2:length(knn)]


#adjust predicted expressions with 10 pcs
#this the function
adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

pxcan_adj <- NULL
for (gene in pxcangenes){
  gene_expr <- en[[gene]]
  pxcan_adj <- cbind(pxcan_adj,adjust_for_covariates(gene_expr, cov_df))
}
pxcan_adj <- as.data.frame(pxcan_adj)
names(pxcan_adj) <- pxcangenes #put back the name of the genes as column

rf_adj <- NULL
for (gene in rfgenes){
  gene_expr <- rf[[gene]]
  rf_adj <- cbind(rf_adj,adjust_for_covariates(gene_expr, cov_df))
}
rf_adj <- as.data.frame(rf_adj)
names(rf_adj) <- rfgenes

svr_adj <- NULL
for (gene in svrgenes){
  gene_expr <- svr[[gene]]
  svr_adj <- cbind(svr_adj, adjust_for_covariates(gene_expr, cov_df))
}
svr_adj <- as.data.frame(svr_adj)
names(svr_adj) <- svrgenes

knn_adj <- NULL
for (gene in knngenes){
  gene_expr <- knn[[gene]]
  knn_adj <- cbind(knn_adj, adjust_for_covariates(gene_expr, cov_df))
}
knn_adj <- as.data.frame(knn_adj)
names(knn_adj) <- knngenes

#put the IID_order back to the adjusted predicted expression
pxcan_adj <- cbind(IID_order, pxcan_adj)
rf_adj <- cbind(IID_order, rf_adj)
svr_adj <- cbind(IID_order, svr_adj)
knn_adj <- cbind(IID_order, knn_adj)

#writeout the adjusted expressions
fwrite(pxcan_adj, file="Z:/data/twas_mesa/en_ALL_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
fwrite(rf_adj, file="Z:/data/twas_mesa/rf_ALL_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
fwrite(svr_adj, file="Z:/data/twas_mesa/svr_ALL_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
fwrite(knn_adj, file="Z:/data/twas_mesa/knn_ALL_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")


#Now run TWAS for each of the pheno

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

#reading the phenotypes with no training samples
#platelet count
plt <- fread(file="Z:/data/twas_mesa/rankplt5_notrain.txt", header=T)
names(plt) <- c("IID","phenotype")

enplt <- inner_join(plt, pxcan_adj, by = c("IID"="IID"))
enplt$IID <- NULL
rfplt <- inner_join(plt, rf_adj, by = c("IID"="IID"))
rfplt$IID <- NULL
svrplt <- inner_join(plt, svr_adj, by = c("IID"="IID"))
svrplt$IID <- NULL
knnplt <- inner_join(plt, knn_adj, by = c("IID"="IID"))
knnplt$IID <- NULL

#do twas
pxcanassoc <- association(enplt, pxcangenes)
#pxcanassoc$p <- as.numeric(pxcanassoc$p)
rfassoc <- association(rfplt, rfgenes)
#rfassoc$p <- as.numeric(rfassoc$p)
svrassoc <- association(svrplt, svrgenes)
#svrassoc$p <- as.numeric(svrassoc$p)
knnassoc <- association(knnplt, knngenes)
#knnassoc$p <- as.numeric(knnassoc$p)

#write out the plt assoc
fwrite(pxcanassoc, file="Z:/data/twas_mesa/en_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")

#read in the assoc again so R can treat the pvalue as numeric
assoc_en <- fread(file="Z:/data/twas_mesa/en_rankplt5_assoc.txt",header=T)
assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rank_hdl_assoc.txt",header=T)
assoc_svr <- fread(file="Z:/data/twas_mesa/svr_rank_hdl_assoc.txt",header=T)
assoc_knn <- fread(file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt",header=T)

#ggplot after filtering
library(qqman)
qq(assoc_en$p, main="Elastic Net Rank_plt5")
qq(assoc_rf$p, main="Random Forest Rank_plt5")
qq(assoc_svr$p, main="Support Vector Rank_plt5")
qq(assoc_knn$p, main="KNN Rank_plt5")
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
#chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
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
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net", annotatePval=0.0005,
          col=c("red", "blue"))

#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest", annotatePval=0.0005,
          col=c("red", "blue"))

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector", annotatePval=0.0005,
          col=c("red", "blue"))

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
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN", annotatePval=0.0005,
          col=c("red", "blue"))




#platelet count
plt <- fread(file="Z:/data/twas_mesa/rankplt5_notrain.txt", header=T)
names(plt) <- c("IID","phenotype")

enplt <- inner_join(plt, pxcan_adj, by = c("IID"="IID"))
enplt$IID <- NULL
rfplt <- inner_join(plt, rf_adj, by = c("IID"="IID"))
rfplt$IID <- NULL
svrplt <- inner_join(plt, svr_adj, by = c("IID"="IID"))
svrplt$IID <- NULL
knnplt <- inner_join(plt, knn_adj, by = c("IID"="IID"))
knnplt$IID <- NULL

#do twas
pxcanassoc <- association(enplt, pxcangenes)
#pxcanassoc$p <- as.numeric(pxcanassoc$p)
rfassoc <- association(rfplt, rfgenes)
#rfassoc$p <- as.numeric(rfassoc$p)
svrassoc <- association(svrplt, svrgenes)
#svrassoc$p <- as.numeric(svrassoc$p)
knnassoc <- association(knnplt, knngenes)
#knnassoc$p <- as.numeric(knnassoc$p)

#write out the plt assoc
fwrite(pxcanassoc, file="Z:/data/twas_mesa/en_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")

#read in the assoc again so R can treat the pvalue as numeric
assoc_en <- fread(file="Z:/data/twas_mesa/en_rankplt5_assoc.txt",header=T)
assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt",header=T)
assoc_svr <- fread(file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt",header=T)
assoc_knn <- fread(file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt",header=T)

#ggplot after filtering
library(qqman)
qq(assoc_en$p, main="Elastic Net Rank_plt5")
qq(assoc_rf$p, main="Random Forest Rank_plt5")
qq(assoc_svr$p, main="Support Vector Rank_plt5")
qq(assoc_knn$p, main="KNN Rank_plt5")
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
#chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
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
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net", annotatePval=0.0005,
          col=c("red", "blue"))

#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest", annotatePval=0.0005,
          col=c("red", "blue"))

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector", annotatePval=0.0005,
          col=c("red", "blue"))

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
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN", annotatePval=0.0005,
          col=c("red", "blue"))





#platelet count
plt <- fread(file="Z:/data/twas_mesa/rankplt5_notrain.txt", header=T)
names(plt) <- c("IID","phenotype")

enplt <- inner_join(plt, pxcan_adj, by = c("IID"="IID"))
enplt$IID <- NULL
rfplt <- inner_join(plt, rf_adj, by = c("IID"="IID"))
rfplt$IID <- NULL
svrplt <- inner_join(plt, svr_adj, by = c("IID"="IID"))
svrplt$IID <- NULL
knnplt <- inner_join(plt, knn_adj, by = c("IID"="IID"))
knnplt$IID <- NULL

#do twas
pxcanassoc <- association(enplt, pxcangenes)
#pxcanassoc$p <- as.numeric(pxcanassoc$p)
rfassoc <- association(rfplt, rfgenes)
#rfassoc$p <- as.numeric(rfassoc$p)
svrassoc <- association(svrplt, svrgenes)
#svrassoc$p <- as.numeric(svrassoc$p)
knnassoc <- association(knnplt, knngenes)
#knnassoc$p <- as.numeric(knnassoc$p)

#write out the plt assoc
fwrite(pxcanassoc, file="Z:/data/twas_mesa/en_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")

#read in the assoc again so R can treat the pvalue as numeric
assoc_en <- fread(file="Z:/data/twas_mesa/en_rankplt5_assoc.txt",header=T)
assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt",header=T)
assoc_svr <- fread(file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt",header=T)
assoc_knn <- fread(file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt",header=T)

#ggplot after filtering
library(qqman)
qq(assoc_en$p, main="Elastic Net Rank_plt5")
qq(assoc_rf$p, main="Random Forest Rank_plt5")
qq(assoc_svr$p, main="Support Vector Rank_plt5")
qq(assoc_knn$p, main="KNN Rank_plt5")
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
#chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
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
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net", annotatePval=0.0005,
          col=c("red", "blue"))

#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest", annotatePval=0.0005,
          col=c("red", "blue"))

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector", annotatePval=0.0005,
          col=c("red", "blue"))

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
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN", annotatePval=0.0005,
          col=c("red", "blue"))





#HDL
plt <- fread(file="Z:/data/twas_mesa/hdl_notrain_rank.txt", header=T)
plt$hdl <- NULL # remove the phentype that is not rank normalize
names(plt) <- c("IID","phenotype")

enplt <- inner_join(plt, pxcan_adj, by = c("IID"="IID"))
enplt$IID <- NULL
rfplt <- inner_join(plt, rf_adj, by = c("IID"="IID"))
rfplt$IID <- NULL
svrplt <- inner_join(plt, svr_adj, by = c("IID"="IID"))
svrplt$IID <- NULL
knnplt <- inner_join(plt, knn_adj, by = c("IID"="IID"))
knnplt$IID <- NULL

#do twas
pxcanassoc <- association(enplt, pxcangenes)
#pxcanassoc$p <- as.numeric(pxcanassoc$p)
rfassoc <- association(rfplt, rfgenes)
#rfassoc$p <- as.numeric(rfassoc$p)
svrassoc <- association(svrplt, svrgenes)
#svrassoc$p <- as.numeric(svrassoc$p)
knnassoc <- association(knnplt, knngenes)
#knnassoc$p <- as.numeric(knnassoc$p)

#write out the plt assoc
fwrite(pxcanassoc, file="Z:/data/twas_mesa/en_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")

#read in the assoc again so R can treat the pvalue as numeric
assoc_en <- fread(file="Z:/data/twas_mesa/en_rankplt5_assoc.txt",header=T)
assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt",header=T)
assoc_svr <- fread(file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt",header=T)
assoc_knn <- fread(file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt",header=T)

#ggplot after filtering
library(qqman)
qq(assoc_en$p, main="Elastic Net Rank_plt5")
qq(assoc_rf$p, main="Random Forest Rank_plt5")
qq(assoc_svr$p, main="Support Vector Rank_plt5")
qq(assoc_knn$p, main="KNN Rank_plt5")
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
#chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
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
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net", annotatePval=0.0005,
          col=c("red", "blue"))

#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest", annotatePval=0.0005,
          col=c("red", "blue"))

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector", annotatePval=0.0005,
          col=c("red", "blue"))

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
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN", annotatePval=0.0005,
          col=c("red", "blue"))
