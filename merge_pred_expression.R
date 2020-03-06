#merge the predicted expression
library(dplyr)
library(tidyr)
library(data.table)

"%&%" = function(a,b) paste (a,b,sep="")

#This is to merge all the ml predicted expressions that was broken in chunks, into one
#Also disregard the variable names, cause I read in and store data in variable regardless of their names

for (pop in c("afa","cau","his")){
  
  print(pop)
  
  
  afa_sam <- fread(file="/home/pokoro/data/lauren_mesa/" %&% pop %&% "_dosages/samples.txt", header=F) # read in sample ID
  names(afa_sam) <- c("FID","IID")
  rf_afa <- NULL
  svr_afa <- NULL
  knn_afa <- NULL
  for (i in 1:22){
    for (j in 1:5){
      rf_afa <- rbind(rf_afa,fread(file="/home/pokoro/data/twas_mesa/ml_pred/" %&% pop %&% "_chr" %&% i %&% "_chunk" %&% j %&% "_rf.txt", header=T))
      svr_afa <- rbind(svr_afa,fread(file="/home/pokoro/data/twas_mesa/ml_pred/" %&% pop %&% "_chr" %&% i %&% "_chunk" %&% j %&% "_svr.txt", header=T))
      knn_afa <- rbind(knn_afa,fread(file="/home/pokoro/data/twas_mesa/ml_pred/" %&% pop %&% "_chr" %&% i %&% "_chunk" %&% j %&% "_knn.txt", header=T))
    }
  }
  
  #remove the gene duplicates
  rf_afa <- subset(rf_afa, !duplicated(rf_afa$gene_id))
  svr_afa <- subset(svr_afa, !duplicated(svr_afa$gene_id))
  knn_afa <- subset(knn_afa, !duplicated(knn_afa$gene_id))
  
  fwrite(rf_afa, file="/home/pokoro/data/twas_mesa/ml_pred/" %&% pop %&% "_allchrom_rf.txt", quote=F,row.names=F,sep="\t")
  fwrite(svr_afa, file="/home/pokoro/data/twas_mesa/ml_pred/" %&% pop %&% "_allchrom_svr.txt", quote=F,row.names=F,sep="\t")
  fwrite(knn_afa, file="/home/pokoro/data/twas_mesa/ml_pred/" %&% pop %&% "_allchrom_knn.txt", quote=F,row.names=F,sep="\t")
  
  #
  #transform the pred exp to the required
  #I used same vairable name for all the data in order to save memory
  # sample_id gene1 gene2 .....
  #rf
  trf <- as.data.frame(t(rf_afa))
  rowrf <- rf_afa$gene_id #take the gene id's
  colnames(trf) <- rowrf
  trf <- trf[-1,] #delete the first row
  trf <- cbind(afa_sam[,2], trf) #merge with the sample id
  fwrite(trf, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_rf.txt", quote=F,row.names=F,sep="\t")
  
  #svr
  trf <- as.data.frame(t(svr_afa))
  rowrf <- rf_afa$gene_id #take the gene id's
  colnames(trf) <- rowrf
  trf <- trf[-1,] #delete the first row
  trf <- cbind(afa_sam[,2], trf) #merge with the sample id
  fwrite(trf, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_svr.txt", quote=F,row.names=F,sep="\t")
  
  #knn
  trf <- as.data.frame(t(knn_afa))
  rowrf <- rf_afa$gene_id #take the gene id's
  colnames(trf) <- rowrf
  trf <- trf[-1,] #delete the first row
  trf <- cbind(afa_sam[,2], trf) #merge with the sample id
  fwrite(trf, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_knn.txt", quote=F,row.names=F,sep="\t")
  
  
}

# merge the afa, cau, and his into one
# 
# for (alg in c("rf","svr","knn")){
#   print(alg)
#   mesa <- NULL
#   for (pop in c("afa","cau","his")){
#     print(pop)
#     mesa <- rbind(mesa, fread(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_" %&% alg %&% ".txt", header=T))
#   }
#   fwrite(mesa, file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_" %&% alg %&% ".txt", quote=F,row.names=F,sep="\t")
# }



#############
#The merging of afa, cau, and his into one as was done above is wrong because, their columns were ordered differently
#Therefore in order to get the correct merge, I need to make the columns to be in the same order
#do this by ordering the cau and his columns to be in same order with afa.
#use the dplyr select method to arrange the columns to be in same order with afa. This is achievable because the columns names
#are exactly the same. The only challenge was that they were in different order.
#Therefore when I read in the cau and his, I will pipe it to the select method, and giving the colnames of afa to it, thereby making it
#to select afa from colnames of cau and his, now because afa cau and his colnames are the same, this will make the cau and his df to be
#ordered the same way afa is ordered. Then we can now rbind all of them into ALL for each of rf, svr, and knn

# afa_rf <- fread(file="Z:/data/twas_mesa/ml_pred/transformed_afa_allchrom_" %&% alg %&% ".txt", header=T, nrows=10) %>% select(colnames(afa_col))
# cau_rf <- fread(file="Z:/data/twas_mesa/ml_pred/transformed_cau_allchrom_" %&% alg %&% ".txt", header=T, nrows=10) %>% select(colnames(afa_col))
# his_rf <- fread(file="Z:/data/twas_mesa/ml_pred/transformed_his_allchrom_" %&% alg %&% ".txt", header=T, nrows=10)
# 
# #so use this loop. first read in afa to take the colnames
# 
# afa_col <- fread(file="Z:/data/twas_mesa/ml_pred/transformed_afa_allchrom_" %&% alg %&% ".txt", header=T, nrows=3)
# 
# for (alg in c("rf","svr","knn")){
#   print(alg)
#   mesa <- NULL
#   for (pop in c("afa","cau","his")){
#     print(pop)
#     afa_col <- fread(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_afa_allchrom_" %&% alg %&% ".txt", header=T, nrows=3)
#     mesa <- rbind(mesa, fread(file="/home/pokoro/data/twas_mesa/ml_pred/transformed_" %&% pop %&% "_allchrom_" %&% alg %&% ".txt", header=T) %>% select(colnames(afa_col)))
#   }
#   fwrite(mesa,file="/home/pokoro/data/twas_mesa/ml_pred/transformed_ALL_allchrom_" %&% alg %&% ".txt", quote=F,row.names=F,sep="\t")
# }




















# #read in the sample IDs
# afa_sam <- fread(file="Z:/data/lauren_mesa/afa_dosages/samples.txt", header=F)
# cau_sam <- fread(file="Z:/data/lauren_mesa/cau_dosages/samples.txt", header=F)
# his_sam <- fread(file="Z:/data/lauren_mesa/his_dosages/samples.txt", header=F)
# 
# #read in the phenotypes and their ID
# trig_mesa_notrain <- fread(file="Z:/data/twas_mesa/trig_notrain.txt", header=T)
# ldl_mesa_notrain <- fread(file="Z:/data/twas_mesa/ldl_train.txt", header=T)
# hdl_mesa_notrain <- fread(file="Z:/data/twas_mesa/hdl_notrain.txt", header=T)
# chol_mesa_notrain <- fread(file="Z:/data/twas_mesa/chol_train.txt", header=T)
# 
##
# merge all the predicted expression chunks and chromosomes per population
# RF
# AFA
# 
# afa_sam <- fread(file="Z:/data/lauren_mesa/afa_dosages/samples.txt", header=F) # read in sample ID
# names(afa_sam) <- c("FID","IID")
# rf_afa <- NULL
# svr_afa <- NULL
# knn_afa <- NULL
# for (i in 1:22){
#   for (j in 1:5){
#     rf_afa <- rbind(rf_afa,fread(file="Z:/data/twas_mesa/ml_pred/afa_chr" %&% i %&% "_chunk" %&% j %&% "_rf.txt", header=T))
#     svr_afa <- rbind(rf_afa,fread(file="Z:/data/twas_mesa/ml_pred/afa_chr" %&% i %&% "_chunk" %&% j %&% "_svr.txt", header=T))
#     knn_afa <- rbind(rf_afa,fread(file="Z:/data/twas_mesa/ml_pred/afa_chr" %&% i %&% "_chunk" %&% j %&% "_knn.txt", header=T))
#   }
# }
# #rf_afa <- fread(file="Z:/data/twas_mesa/ml_pred/afa_chr21_chunk1_svr.txt", header=T)
# fwrite(rf_afa, file="Z:/data/twas_mesa/ml_pred/afa_allchrom_rf.txt", quote=F,row.names=F,sep="\t")
# fwrite(svr_afa, file="Z:/data/twas_mesa/ml_pred/afa_allchrom_svr.txt", quote=F,row.names=F,sep="\t")
# fwrite(knn_afa, file="Z:/data/twas_mesa/ml_pred/afa_allchrom_knn.txt", quote=F,row.names=F,sep="\t")
# 
# #
# #transform the pred exp to the required
# # sample_id gene1 gene2 .....
# #rf
# trf <- as.data.frame(t(rf_afa))
# rowrf <- rf_afa$gene_id #take the gene id's
# colnames(trf) <- rowrf
# trf <- trf[-1,] #delete the first row
# trf <- cbind(afa_sam[,2], trf) #merge with the sample id
# fwrite(trf, file="Z:/data/twas_mesa/ml_pred/transformed_afa_allchrom_rf.txt", quote=F,row.names=F,sep="\t")
# 
# #svr
# trf <- as.data.frame(t(svr_afa))
# rowrf <- rf_afa$gene_id #take the gene id's
# colnames(trf) <- rowrf
# trf <- trf[-1,] #delete the first row
# trf <- cbind(afa_sam[,2], trf) #merge with the sample id
# fwrite(trf, file="Z:/data/twas_mesa/ml_pred/transformed_afa_allchrom_svr.txt", quote=F,row.names=F,sep="\t")
# 
# #knn
# trf <- as.data.frame(t(knn_afa))
# rowrf <- rf_afa$gene_id #take the gene id's
# colnames(trf) <- rowrf
# trf <- trf[-1,] #delete the first row
# trf <- cbind(afa_sam[,2], trf) #merge with the sample id
# fwrite(trf, file="Z:/data/twas_mesa/ml_pred/transformed_afa_allchrom_knn.txt", quote=F,row.names=F,sep="\t")
