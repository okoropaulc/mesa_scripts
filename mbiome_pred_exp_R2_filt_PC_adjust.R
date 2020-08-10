library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
library(dplyr)
library(tidyr)

#loop through
######################################################################
#read in the best grids results for the Pops


all_rfgrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
all_svrgrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
all_knngrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)

#read in the predicted expression
all_pred_rf <- fread(file="Z:/data/mesa_models/python_ml_models/results/ALL_2_METS_rf_predicted_expr_allchrom.txt", header=T)
all_pred_svr <- fread(file="Z:/data/mesa_models/python_ml_models/results/ALL_2_METS_svr_rbf_predicted_expr_allchrom.txt", header=T)
all_pred_knn <- fread(file="Z:/data/mesa_models/python_ml_models/results/ALL_2_METS_knn_predicted_expr_allchrom.txt", header=T)


#convert to gene_id and cv_r2 to characters
all_rfgrid$Gene_ID <- as.character(all_rfgrid$Gene_ID)
all_rfgrid$CV_R2 <- as.numeric(all_rfgrid$CV_R2)
all_svrgrid$Gene_ID <- as.character(all_svrgrid$Gene_ID)
all_svrgrid$CV_R2 <- as.character(all_svrgrid$CV_R2)
all_knngrid$Gene_ID <- as.character(all_knngrid$Gene_ID)
all_knngrid$CV_R2 <- as.numeric(all_knngrid$CV_R2)

#first remove decimals from the gene id. the number of genes=9623 is same for all algs
for (i in 1:nrow(all_rfgrid)){
  all_rfgrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_rfgrid$Gene_ID[i])
} #just to remove the decimal places in the gene_id

for (i in 1:nrow(all_svrgrid)){
  all_svrgrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_svrgrid$Gene_ID[i])
}

for (i in 1:nrow(all_knngrid)){
  all_knngrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_knngrid$Gene_ID[i])
}

#filter by cv R2 > 0.01
all_rfgrid <- subset(all_rfgrid, CV_R2 > 0.01)
all_rfgrid <- all_rfgrid[,c(1,2,3)]# keep only gene_id, gene_name, and cv_r2
all_svrgrid <- subset(all_svrgrid, CV_R2 > 0.01)
all_svrgrid <- all_svrgrid[,c(1,2,3)]
all_knngrid <- subset(all_knngrid, CV_R2 > 0.01)
all_knngrid <- all_knngrid[,c(1,2,3)]

#keep predicted genes that have cv_r2 > 0.01
IID <- all_pred_rf$V1
all_pred_rf <- all_pred_rf[,-1]
rfcol <- colnames(all_pred_rf)
trf_pred <- as.data.frame(t(all_pred_rf))
names(trf_pred) <- IID
trf_pred <- cbind.data.frame(rfcol, trf_pred)
trf_pred <- inner_join(all_rfgrid, trf_pred, by = c("Gene_ID"="rfcol"))
trf_pred <- trf_pred[,-c(2,3)]
rem_genes <- trf_pred$Gene_ID
rf <- as.data.frame(t(trf_pred))
rf <- rf[-1,]
names(rf) <- rem_genes
rf <- cbind.data.frame(IID, rf)

#rf <- all_pred_rf %>% select(c(trf_pred$Gene_ID))


IID <- all_pred_svr$V1
all_pred_svr <- all_pred_svr[,-1]
svrcol <- colnames(all_pred_svr)
tsvr_pred <- as.data.frame(t(all_pred_svr))
names(tsvr_pred) <- IID
tsvr_pred <- cbind(svrcol, tsvr_pred)
tsvr_pred <- inner_join(all_svrgrid, tsvr_pred, by = c("Gene_ID"="svrcol"))
tsvr_pred <- tsvr_pred[,-c(2,3)]
rem_genes <- tsvr_pred$Gene_ID
svr <- as.data.frame(t(tsvr_pred))
svr <- svr[-1,]
names(svr) <- rem_genes
svr <- cbind(IID, svr)

#svr <- all_pred_svr %>% select(c("V1", all_svrgrid$Gene_ID))


IID <- all_pred_knn$V1
all_pred_knn <- all_pred_knn[,-1]
knncol <- colnames(all_pred_knn)
tknn_pred <- as.data.frame(t(all_pred_knn))
names(tknn_pred) <- IID
tknn_pred <- cbind(knncol, tknn_pred)
tknn_pred <- inner_join(all_knngrid, tknn_pred, by = c("Gene_ID"="knncol"))
tknn_pred <- tknn_pred[,-c(2,3)]
rem_genes <- tknn_pred$Gene_ID
knn <- as.data.frame(t(tknn_pred))
knn <- knn[-1,]
names(knn) <- rem_genes
knn <- cbind(IID, knn)

#knn <- all_pred_knn %>% select(c("V1", all_knngrid$Gene_ID))

#write out the filtered gene expression
#fwrite(rf, file="Z:/paul/mbiome_assoc/ALL_2_METS_rf_pred_exp_filteredbyR2_allchrom.txt", row.names=F, quote=F, sep="\t")
#fwrite(svr, file="Z:/paul/mbiome_assoc/ALL_2_METS_svr_pred_exp_filteredbyR2_allchrom.txt", row.names=F, quote=F, sep="\t")
#fwrite(knn, file="Z:/paul/mbiome_assoc/ALL_2_METS_knn_pred_exp_filteredbyR2_allchrom.txt", row.names=F, quote=F, sep="\t")

#read in again to stop the numerics from being factor
rf <- fread(file="Z:/paul/mbiome_assoc/ALL_2_METS_rf_pred_exp_filteredbyR2_allchrom.txt", header=T)
svr <- fread(file="Z:/paul/mbiome_assoc/ALL_2_METS_svr_pred_exp_filteredbyR2_allchrom.txt", header=T)
knn <- fread(file="Z:/paul/mbiome_assoc/ALL_2_METS_knn_pred_exp_filteredbyR2_allchrom.txt", header=T)



# #do same filtering for EN. use en_all as varibale name


en_all <- fread(file="Z:/data/mesa_models/split_mesa/results/all_chr_ALL_model_summaries.txt", header=T)

en_all$gene_id <- as.character(en_all$gene_id)
for (i in 1:nrow(en_all)){
  en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
} #just to remove the decimal places in the gene_id


en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)

#filter by cv R2 > 0.01
en_all <- subset(en_all, cv_R2_avg > 0.01)
en_all <- en_all[,c(1,2,10)]# keep only gene_id, gene_name, and cv_r2


#read in the EN ALL predicted expression

all_pred_en <- fread(file = "Z:/data/mesa_models/mets_dosages/ALL_2_METS_predicted_expression.txt", header = T)

#remove the decimal in the gene id
for (i in 3:ncol(all_pred_en)){
  colnames(all_pred_en)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_en)[i])
}

#keep predicted genes that have cv_r2 >= 0.01
en <- all_pred_en %>% select(c("IID", en_all$gene_id))

#write out
fwrite(en, file="Z:/paul/mbiome_assoc/ALL_2_METS_en_pred_exp_filteredbyR2_allchrom.txt", row.names=F, quote=F, sep="\t")
#read in again
en <- fread(file="Z:/paul/mbiome_assoc/ALL_2_METS_en_pred_exp_filteredbyR2_allchrom.txt", header=T)




##Adjust by PC

adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}


#read in the pcs of the ld pruned mesa ALL
mets_pc <- fread(file="Z:/data/METS_model/mets_cov_pca_10_eigenvec.txt", header=T)

for (i in 1:nrow(mets_pc)){
  mets_pc$id[i] <- paste(mets_pc$FID[i], mets_pc$IID[i], sep="_") 
}

#Now arrange the sample IID of the predicted expressions to be in same order with the mesa_pc sample IID. arrange pheno with this too
#use inner_join
IID_order <- data.frame(IID=mets_pc[,"id"])
names(IID_order) <- ("IID")
#fwrite(IID_order, file="/home/pokoro/data/twas_mesa/samples_for_adj_pred_exp.txt", row.names=F, quote=F)

rf <- inner_join(IID_order, rf, by = c("IID"="IID"))
svr <- inner_join(IID_order, svr, by = c("IID"="IID"))
knn <- inner_join(IID_order, knn, by = c("IID"="IID"))
en <- inner_join(IID_order, en, by = c("IID"="IID"))

cov_df <- mets_pc[,3:12]#take only the 10 pcs since the pc sample IID is in same order with the predicted expression

#take the genes names in the predicted expression dataframe, so I can use it to iteratively adjust the expression with pc
pxcangenes <- colnames(en)[2:length(en)]
rfgenes <- colnames(rf)[2:length(rf)]
svrgenes <- colnames(svr)[2:length(svr)]
knngenes <- colnames(knn)[2:length(knn)]



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
#fwrite(pxcan_adj, file="/home/pokoro/data/twas_mesa/en_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
#fwrite(rf_adj, file="/home/pokoro/data/twas_mesa/rf_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
#fwrite(svr_adj, file="/home/pokoro/data/twas_mesa/svr_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
#fwrite(knn_adj, file="/home/pokoro/data/twas_mesa/knn_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")




#
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
#read in the mbiome diversity index file
diversity <- read.table(file="Z:/mets61/mets61_pheno_alpha_D_index_sorted.txt", header=T)
diversity$Obesity_status <- as.character(diversity$Obesity_status)
#diversity$obese <- ""

for (i in 1:nrow(diversity)){
  #  if(diversity$Obesity_status[i]=="Lean"){
  #    diversity$obese[i] <- 0
  #  }
  ifelse(diversity$Obesity_status[i]=="Lean", diversity$obese[i]<-0, diversity$obese[i]<-1)
}

#diversity index shannon
shannon <- diversity %>% select("V2","Shannon_index")
names(shannon) <- c("IID","phenotype")

enshannon <- inner_join(shannon, pxcan_adj, by = c("IID"="IID"))
enshannon$IID <- NULL
rfshannon <- inner_join(shannon, rf_adj, by = c("IID"="IID"))
rfshannon$IID <- NULL
svrshannon <- inner_join(shannon, svr_adj, by = c("IID"="IID"))
svrshannon$IID <- NULL
knnshannon <- inner_join(shannon, knn_adj, by = c("IID"="IID"))
knnshannon$IID <- NULL

#do twas
pxcanassoc <- association(enshannon, pxcangenes)
rfassoc <- association(rfshannon, rfgenes)
svrassoc <- association(svrshannon, svrgenes)
knnassoc <- association(knnshannon, knngenes)


#writeout the assoc results
fwrite(pxcanassoc, file="Z:/paul/mbiome_assoc/en_mets_shannon_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/paul/mbiome_assoc/rf_mets_shannon_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/paul/mbiome_assoc/svr_mets_shannon_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/paul/mbiome_assoc/knn_mets_shannon_twas.txt", row.names=F, quote=F, sep="\t")

#check 
pxcanassoc <- fread(file="Z:/paul/mbiome_assoc/en_mets_shannon_twas.txt", header=T)
rfassoc <- fread(file="Z:/paul/mbiome_assoc/rf_mets_shannon_twas.txt", header=T)
svrassoc <- fread(file="Z:/paul/mbiome_assoc/svr_mets_shannon_twas.txt", header=T)
knnassoc <- fread(file="Z:/paul/mbiome_assoc/knn_mets_shannon_twas.txt", header=T)


#diversity index Simpson
simpson <- diversity %>% select("V2","Simpson_index")
names(simpson) <- c("IID","phenotype")

ensimpson <- inner_join(simpson, pxcan_adj, by = c("IID"="IID"))
ensimpson$IID <- NULL
rfsimpson <- inner_join(simpson, rf_adj, by = c("IID"="IID"))
rfsimpson$IID <- NULL
svrsimpson <- inner_join(simpson, svr_adj, by = c("IID"="IID"))
svrsimpson$IID <- NULL
knnsimpson <- inner_join(simpson, knn_adj, by = c("IID"="IID"))
knnsimpson$IID <- NULL

#do twas
pxcanassoc <- association(ensimpson, pxcangenes)
rfassoc <- association(rfsimpson, rfgenes)
svrassoc <- association(svrsimpson, svrgenes)
knnassoc <- association(knnsimpson, knngenes)


#writeout the assoc results
fwrite(pxcanassoc, file="Z:/paul/mbiome_assoc/en_mets_simpson_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/paul/mbiome_assoc/rf_mets_simpson_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/paul/mbiome_assoc/svr_mets_simpson_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/paul/mbiome_assoc/knn_mets_simpson_twas.txt", row.names=F, quote=F, sep="\t")

#check 
pxcanassoc <- fread(file="Z:/paul/mbiome_assoc/en_mets_simpson_twas.txt", header=T)
rfassoc <- fread(file="Z:/paul/mbiome_assoc/rf_mets_simpson_twas.txt", header=T)
svrassoc <- fread(file="Z:/paul/mbiome_assoc/svr_mets_simpson_twas.txt", header=T)
knnassoc <- fread(file="Z:/paul/mbiome_assoc/knn_mets_simpson_twas.txt", header=T)




#obesity obese=1 lean=0
obese <- diversity %>% select("V2","obese")
names(obese) <- c("IID","phenotype")

enobese <- inner_join(obese, pxcan_adj, by = c("IID"="IID"))
enobese$IID <- NULL
rfobese <- inner_join(obese, rf_adj, by = c("IID"="IID"))
rfobese$IID <- NULL
svrobese <- inner_join(obese, svr_adj, by = c("IID"="IID"))
svrobese$IID <- NULL
knnobese <- inner_join(obese, knn_adj, by = c("IID"="IID"))
knnobese$IID <- NULL

#do twas
pxcanassoc <- association(enobese, pxcangenes, test_type = "logistic")
rfassoc <- association(rfobese, rfgenes, test_type = "logistic")
svrassoc <- association(svrobese, svrgenes, test_type = "logistic")
knnassoc <- association(knnobese, knngenes, test_type = "logistic")


#writeout the assoc results
fwrite(pxcanassoc, file="Z:/paul/mbiome_assoc/en_mets_obese_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/paul/mbiome_assoc/rf_mets_obese_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/paul/mbiome_assoc/svr_mets_obese_twas.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/paul/mbiome_assoc/knn_mets_obese_twas.txt", row.names=F, quote=F, sep="\t")

#check 
pxcanassoc <- fread(file="Z:/paul/mbiome_assoc/en_mets_obese_twas.txt", header=T)
rfassoc <- fread(file="Z:/paul/mbiome_assoc/rf_mets_obese_twas.txt", header=T)
svrassoc <- fread(file="Z:/paul/mbiome_assoc/svr_mets_obese_twas.txt", header=T)
knnassoc <- fread(file="Z:/paul/mbiome_assoc/knn_mets_obese_twas.txt", header=T)



######
#Further Adjust the PC adjusted predicted expression with diversity index
#Shannon
shannon <- diversity %>% select("V2","Shannon_index")
names(shannon) <- c("IID","shannon")
shannon <- inner_join(IID_order, shannon, by = c("IID"="IID"))

pxcan_adj <- inner_join(shannon, pxcan_adj, by = c("IID"="IID"))
rf_adj <- inner_join(shannon, rf_adj, by = c("IID"="IID"))
svr_adj <- inner_join(shannon, svr_adj, by = c("IID"="IID"))
knn_adj <- inner_join(shannon, knn_adj, by = c("IID"="IID"))

metsIID <- rf_adj$IID #Just to keep the IID order

shannon$IID <- NULL

pxcan_shannon_adj <- NULL
for (gene in pxcangenes){
  gene_expr <- pxcan_adj[[gene]]
  pxcan_shannon_adj <- cbind(pxcan_shannon_adj,adjust_for_covariates(gene_expr, shannon))
}

pxcan_shannon_adj <- as.data.frame(pxcan_shannon_adj)
names(pxcan_shannon_adj) <- pxcangenes #put back the name of the genes as column

rf_shannon_adj <- NULL
for (gene in rfgenes){
  gene_expr <- rf_adj[[gene]]
  rf_shannon_adj <- cbind(rf_shannon_adj,adjust_for_covariates(gene_expr, shannon))
}

rf_shannon_adj <- as.data.frame(rf_shannon_adj)
names(rf_shannon_adj) <- rfgenes

svr_shannon_adj <- NULL
for (gene in svrgenes){
  gene_expr <- svr_adj[[gene]]
  svr_shannon_adj <- cbind(svr_shannon_adj, adjust_for_covariates(gene_expr, shannon))
}

svr_shannon_adj <- as.data.frame(svr_shannon_adj)
names(svr_shannon_adj) <- svrgenes

knn_shannon_adj <- NULL
for (gene in knngenes){
  gene_expr <- knn_adj[[gene]]
  knn_shannon_adj <- cbind(knn_shannon_adj, adjust_for_covariates(gene_expr, shannon))
}

knn_shannon_adj <- as.data.frame(knn_shannon_adj)
names(knn_shannon_adj) <- knngenes

#put the IID_order back to the shannon_adjusted predicted expression
pxcan_shannon_adj <- cbind(metsIID, pxcan_shannon_adj)
rf_shannon_adj <- cbind(metsIID, rf_shannon_adj)
svr_shannon_adj <- cbind(metsIID, svr_shannon_adj)
knn_shannon_adj <- cbind(metsIID, knn_shannon_adj)



####
#Now do assoication for Obesity

#obesity obese=1 lean=0
obese <- diversity %>% select("V2","obese")
names(obese) <- c("IID","phenotype")

enobese <- inner_join(obese, pxcan_shannon_adj, by = c("IID"="metsIID"))
enobese$IID <- NULL
rfobese <- inner_join(obese, rf_shannon_adj, by = c("IID"="metsIID"))
rfobese$IID <- NULL
svrobese <- inner_join(obese, svr_shannon_adj, by = c("IID"="metsIID"))
svrobese$IID <- NULL
knnobese <- inner_join(obese, knn_shannon_adj, by = c("IID"="metsIID"))
knnobese$IID <- NULL

#do twas
pxcanassoc <- association(enobese, pxcangenes, test_type = "logistic")
rfassoc <- association(rfobese, rfgenes, test_type = "logistic")
svrassoc <- association(svrobese, svrgenes, test_type = "logistic")
knnassoc <- association(knnobese, knngenes, test_type = "logistic")


#writeout the assoc results
fwrite(pxcanassoc, file="Z:/paul/mbiome_assoc/en_mets_obese_twas_adjusted_with_shannon.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/paul/mbiome_assoc/rf_mets_obese_twas_adjusted_with_shannon.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/paul/mbiome_assoc/svr_mets_obese_twas_adjusted_with_shannon.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/paul/mbiome_assoc/knn_mets_obese_twas_adjusted_with_shannon.txt", row.names=F, quote=F, sep="\t")

#check 
pxcanassoc <- fread(file="Z:/paul/mbiome_assoc/en_mets_obese_twas_adjusted_with_shannon.txt", header=T)
rfassoc <- fread(file="Z:/paul/mbiome_assoc/rf_mets_obese_twas_adjusted_with_shannon.txt", header=T)
svrassoc <- fread(file="Z:/paul/mbiome_assoc/svr_mets_obese_twas_adjusted_with_shannon.txt", header=T)
knnassoc <- fread(file="Z:/paul/mbiome_assoc/knn_mets_obese_twas_adjusted_with_shannon.txt", header=T)
