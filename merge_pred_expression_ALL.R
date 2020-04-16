library(dplyr)
library(tidyr)
library(data.table)

"%&%" = function(a,b) paste (a,b,sep="")

#############
#The merging of afa, cau, and his into one as was done in previous script is wrong because, their columns were ordered differently
#Therefore in order to get the correct merge, I need to make the columns to be in the same order
#do this by ordering the cau and his columns to be in same order with afa.
#use the dplyr select method to arrange the columns to be in same order with afa. This is achievable because the columns names
#are exactly the same. The only challenge was that they were in different order.
#Therefore when I read in the cau and his, I will pipe it to the select method, and giving the colnames of afa to it, thereby making it
#to select afa from colnames of cau and his, now because afa cau and his colnames are the same, this will make the cau and his df to be
#ordered the same way afa is ordered. Then we can now rbind all of them into ALL for each of rf, svr, and knn

##so use this loop. first read in afa to take the colnames

for (alg in c("rf","svr","knn")){
  print(alg)
  mesa <- NULL
  for (pop in c("afa","cau","his")){
    print(pop)
    afa_col <- fread(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_afa_allchrom_" %&% alg %&% ".txt", header=T, nrows=3)
    mesa <- rbind(mesa, fread(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_" %&% alg %&% ".txt", header=T) %>% select(colnames(afa_col)))
  }
  fwrite(mesa,file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_" %&% alg %&% ".txt", quote=F,row.names=F,sep="\t")
}


#merge the pops of en
# afa_en <- fread(file="Z:/data/twas_mesa/AFA_en_predicted_expression.txt", header = T, nrows=10)
# cau_en <- fread(file="Z:/data/twas_mesa/CAU_en_predicted_expression.txt", header = T, nrows=10)
# his_en <- fread(file="Z:/data/twas_mesa/HIS_en_predicted_expression.txt", header = T, nrows=10)

#loop it and save memory
#read in afa cols and use it to arrange the other pops col
# mesa <- NULL
# col_ord <- fread(file="/home/pokoro/data/twas_mesa/AFA_en_predicted_expression.txt", header = T, nrows=3)
# for (i in c("AFA", "CAU", "HIS")){
#   print(i)
#   mesa <- rbind(mesa, fread(file="/home/pokoro/data/twas_mesa/" %&% i %&% "_en_predicted_expression.txt", header = T) %>% select(colnames(col_ord)))
# }
# fwrite(mesa,file="/home/pokoro/data/twas_mesa/ALL_en_predicted_expression.txt", quote=F,row.names=F,sep="\t")










######################################################################
#read in the best gird results for the ALL
all_rfgrid <- fread(file="/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
all_svrgrid <- fread(file="/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
all_knngrid <- fread(file="/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)

#read in the predicted expression
all_pred_rf <- fread(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_rf.txt", header=T)
all_pred_svr <- read.table(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_svr.txt", header=T)
all_pred_knn <- read.table(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_knn.txt", header=T)

#remove the decimals in the gene id (colnames) of the predicted expression
for (i in 2:9624){
  colnames(all_pred_rf)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_rf)[i])
  colnames(all_pred_svr)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_svr)[i])
  colnames(all_pred_knn)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_knn)[i])
}

#convert to gene_id and cv_r2 to characters
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
all_rfgrid <- subset(all_rfgrid, CV_R2 >= 0.01)
all_rfgrid <- all_rfgrid[,c(1,2,3)]# keep only gene_id, gene_name, and cv_r2
all_svrgrid <- subset(all_svrgrid, CV_R2 >= 0.01)
all_svrgrid <- all_svrgrid[,c(1,2,3)]
all_knngrid <- subset(all_knngrid, CV_R2 >= 0.01)
all_knngrid <- all_knngrid[,c(1,2,3)]

#keep predicted genes that have cv_r2 >= 0.01
rf <- all_pred_rf %>% select(c("IID", all_rfgrid$Gene_ID))
svr <- all_pred_svr %>% select(c("IID", all_svrgrid$Gene_ID))
knn <- all_pred_knn %>% select(c("IID", all_knngrid$Gene_ID))

#write out the filtered gene expression
fwrite(rf, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_rf.txt", row.names=F, quote=F, sep="\t")
fwrite(svr, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_svr.txt", row.names=F, quote=F, sep="\t")
fwrite(knn, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_knn.txt", row.names=F, quote=F, sep="\t")



# #do same filtering for EN
# #file does not have header, so read in ones with header and copy the header
# en_cau <- fread(file="/home/pokoro/data/mesa_models/MESA.CAU.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=T)
# en_all <- fread(file="/home/pokoro/no_header_MESA.ALL.WG.PC3.PF10.1000kb_WG_model_summaries.txt", header=F)
# header <- colnames(en_cau)
# names(en_all) <- header
# 
# en_all$gene_id <- as.character(en_all$gene_id)
for (i in 1:nrow(en_all)){
  en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
} #just to remove the decimal places in the gene_id
# 
# 
# en_all$cv_R2_avg <- as.numeric(en_all$cv_R2_avg)
# 
# #filter by cv R2 > 0.01
# en_all <- subset(en_all, cv_R2_avg >= 0.01)
# en_all <- en_all[,c(1,2,10)]# keep only gene_id, gene_name, and cv_r2
# 
# 
# #read in the EN ALL predicted expression
# 
# all_pred_en <- fread(file="/home/pokoro/data/twas_mesa/ALL_en_predicted_expression.txt", header=T)
# 
# #remove the decimal in the gene id
for (i in 3:7808){
  colnames(all_pred_en)[i] <- gsub('\\.[0-9]+','',colnames(all_pred_en)[i])
}
# 
# #keep predicted genes that have cv_r2 >= 0.01
# en <- all_pred_en %>% select(c("IID", en_all$gene_id))
# 
# #write out
# fwrite(en, file="/home/pokoro/data/twas_mesa/ALL_en_predicted_expression_r2_0.01.txt", row.names=F, quote=F, sep="\t")


# #check filtered predicted expression
# rf <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_rf.txt", header = T, nrows = 10)
# svr <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_svr.txt", header = T, nrows = 10)
# knn <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_ALL_allchrom_r2_0.01_knn.txt", header = T, nrows = 10)
# en <- fread(file = "Z:/data/twas_mesa/ALL_en_predicted_expression_r2_0.01.txt", header = T, nrows = 10)
