library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
library(dplyr)
library(tidyr)

#loop through
######################################################################
#read in the best grids results for the Pops

# for (pop in c("afa", "cau", "his")){
#   
#   all_rfgrid <- fread(file="/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/" %&% toupper(pop) %&% "_best_grid_rf_all_chrom.txt", header=T)
#   all_svrgrid <- fread(file="/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/" %&% toupper(pop) %&% "_best_grid_svr_all_chrom.txt", header=T)
#   all_knngrid <- fread(file="/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/" %&% toupper(pop) %&% "_best_grid_knn_all_chrom.txt", header=T)
#   
#   #read in the predicted expression
#   all_pred_rf <- fread(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_rf.txt", header=T)
#   all_pred_svr <- read.table(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_svr.txt", header=T)
#   all_pred_knn <- read.table(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_knn.txt", header=T)
#   
#   #remove the decimals in the gene id (colnames) of the predicted expression
#   for (i in 2:ncol(all_pred_rf)){
#     colnames(all_pred_rf)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_rf)[i])
#   }
#   
#   for (i in 2:ncol(all_pred_svr)){
#     colnames(all_pred_svr)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_svr)[i])
#   }
#   
#   for (i in 2:ncol(all_pred_knn)){
#     colnames(all_pred_knn)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_knn)[i])
#   }
#   
#   #convert to gene_id and cv_r2 to characters
#   all_rfgrid$Gene_ID <- as.character(all_rfgrid$Gene_ID)
#   all_rfgrid$CV_R2 <- as.numeric(all_rfgrid$CV_R2)
#   all_svrgrid$Gene_ID <- as.character(all_svrgrid$Gene_ID)
#   all_svrgrid$CV_R2 <- as.character(all_svrgrid$CV_R2)
#   all_knngrid$Gene_ID <- as.character(all_knngrid$Gene_ID)
#   all_knngrid$CV_R2 <- as.numeric(all_knngrid$CV_R2)
#   
#   #first remove decimals from the gene id. the number of genes=9623 is same for all algs
#   for (i in 1:nrow(all_rfgrid)){
#     all_rfgrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_rfgrid$Gene_ID[i])
#   } #just to remove the decimal places in the gene_id
#   
#   for (i in 1:nrow(all_svrgrid)){
#     all_svrgrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_svrgrid$Gene_ID[i])
#   }
#   
#   for (i in 1:nrow(all_knngrid)){
#     all_knngrid$Gene_ID[i] <- gsub('\\.[0-9]+','',all_knngrid$Gene_ID[i])
#   }
#   #filter by cv R2 > 0.01
#   all_rfgrid <- subset(all_rfgrid, CV_R2 > 0.01)
#   all_rfgrid <- all_rfgrid[,c(1,2,3)]# keep only gene_id, gene_name, and cv_r2
#   all_svrgrid <- subset(all_svrgrid, CV_R2 > 0.01)
#   all_svrgrid <- all_svrgrid[,c(1,2,3)]
#   all_knngrid <- subset(all_knngrid, CV_R2 > 0.01)
#   all_knngrid <- all_knngrid[,c(1,2,3)]
#   
#   #keep predicted genes that have cv_r2 >= 0.01
#   rf <- all_pred_rf %>% select(c("IID", all_rfgrid$Gene_ID))
#   svr <- all_pred_svr %>% select(c("IID", all_svrgrid$Gene_ID))
#   knn <- all_pred_knn %>% select(c("IID", all_knngrid$Gene_ID))
#   
#   #write out the filtered gene expression
#   fwrite(rf, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_r2_0.01_rf.txt", row.names=F, quote=F, sep="\t")
#   fwrite(svr, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_r2_0.01_svr.txt", row.names=F, quote=F, sep="\t")
#   fwrite(knn, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_r2_0.01_knn.txt", row.names=F, quote=F, sep="\t")
#   
# }
# 
# 
# 
# # #do same filtering for EN. use en_all as varibale name
# 
# for (pop in c("AFA", "CAU", "HIS")){
#   
#   en_all <- fread(file="/home/pokoro/data/mesa_models/split_mesa/results/all_chr_" %&% pop %&% "_model_summaries.txt", header=T)
#   
#   en_all$gene_id <- as.character(en_all$gene_id)
#   for (i in 1:nrow(en_all)){
#     en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
#   } #just to remove the decimal places in the gene_id
#   
#   
#   en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)
#   
#   #filter by cv R2 > 0.01
#   en_all <- subset(en_all, cv_R2_avg > 0.01)
#   en_all <- en_all[,c(1,2,10)]# keep only gene_id, gene_name, and cv_r2
#   
#   
#   #read in the EN ALL predicted expression
#   
#   all_pred_en <- fread(file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_en_predicted_expression.txt", header=T)
#   
#   #remove the decimal in the gene id
#   for (i in 3:ncol(all_pred_en)){
#     colnames(all_pred_en)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_en)[i])
#   }
#   
#   #keep predicted genes that have cv_r2 >= 0.01
#   en <- all_pred_en %>% select(c("IID", en_all$gene_id))
#   
#   #write out
#   fwrite(en, file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_en_predicted_expression_r2_0.01.txt", row.names=F, quote=F, sep="\t")
#   
#   
# }
# 





###### Use this to filter for predixcan predicted expression for AFA CAU HIS

for (pop in c("AFA", "CAU", "HIS")){
  
  en_all <- fread(file="/home/pokoro/data/mesa_models/split_mesa/results/all_chr_" %&% pop %&% "_model_summaries.txt", header=T)
  en_all$gene_id <- as.character(en_all$gene_id)
  for (i in 1:nrow(en_all)){
    en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
  } #
  en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)
  
  #filter by cv R2 > 0.01
  en_all <- subset(en_all, cv_R2_avg > 0.01)
  en_all <- en_all[,c(1,2,10)]# keep only gene_id, gene_name, and cv_r2
  
  
  #read in the EN ALL predicted expression
  
  all_pred_en <- fread(file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_en_predicted_expression.txt", header=T)
  
  #remove the decimal in the gene id
  for (i in 3:ncol(all_pred_en)){
    colnames(all_pred_en)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_en)[i])
  }
  
  #keep predicted genes that have cv_r2 >= 0.01
  ten <- (all_pred_en)
  ten$FID <- NULL
  IID <- ten$IID
  ten$IID <- NULL
  ten <- as.data.frame(t(ten))
  genes <- rownames(ten)
  ten <- cbind(genes,ten)
  ten$genes <- as.character(ten$genes)
  en <- inner_join(en_all, ten, by = c("gene_id"="genes"))
  r2genes <- en$gene_id
  en[,c(1,2,3)] <- NULL
  en <- as.data.frame(t(en))
  names(en) <- r2genes
  en <- cbind(IID,en)
  #en <- all_pred_en %>% select(c("IID", en_all$gene_id))
  
  #write out
  fwrite(en, file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_en_predicted_expression_r2_0.01.txt", row.names=F, quote=F, sep="\t")
  
  
}

# #check filtered predicted expression
# rf <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_rf.txt", header = T, nrows = 10)
# svr <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_svr.txt", header = T, nrows = 10)
# knn <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_knn.txt", header = T, nrows = 10)
# en <- fread(file = "Z:/data/twas_mesa/ALL_en_predicted_expression_r2_0.01.txt", header = T, nrows = 10)
