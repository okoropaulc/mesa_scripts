#This script will use the pcs to adjust the predicted expression
#it will also run twas of the expression and phenotype

library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
library(dplyr)


for (pop in c("AFA", "CAU", "HIS")){
  
  #read in the pcs of the ld pruned mesa ALL
  mesa_pc <- read.table(file="/home/pokoro/data/lauren_mesa/ALL/" %&% tolower(pop) %&% "/" %&% tolower(pop) %&% "_pca.eigenvec", header=T)
  
  #read in the predicted expression already filtered by r2 of each algorithm
  rf <- fread(file = "/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% tolower(pop) %&% "_allchrom_r2_0.01_rf.txt", header = T)
  svr <- fread(file = "/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% tolower(pop) %&% "_allchrom_r2_0.01_svr.txt", header = T)
  knn <- fread(file = "/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% tolower(pop) %&% "_allchrom_r2_0.01_knn.txt", header = T)
  en <- fread(file = "/home/pokoro/data/twas_mesa/" %&% pop %&% "_en_predicted_expression_r2_0.01.txt", header = T)
  
  #keep copy of the files, just to keep a copy of these files because it takes long to read them in rstudio
  rf_cp <-rf
  svr_cp <- svr
  knn_cp <- knn
  en_cp <- en
  
  #Now arrange the sample IID of the predicted expressions to be in same order with the mesa_pc sample IID. arrange pheno with this too
  #use inner_join
  IID_order <- data.frame(IID=mesa_pc[,2])
  #fwrite(IID_order, file="/home/pokoro/data/twas_mesa/samples_for_adj_pred_exp.txt", row.names=F, quote=F)
  
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
  #this is the function
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
  fwrite(pxcan_adj, file="/home/pokoro/data/twas_mesa/en_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
  fwrite(rf_adj, file="/home/pokoro/data/twas_mesa/rf_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
  fwrite(svr_adj, file="/home/pokoro/data/twas_mesa/svr_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
  fwrite(knn_adj, file="/home/pokoro/data/twas_mesa/knn_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", row.names=F, quote=F, sep="\t")
  
  
}
