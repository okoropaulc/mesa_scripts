#DO the spearman correlation of the mesa2mets predicted expression vs mets measured expression
# Filter by cvR2>0.01
library(data.table)
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)

#read in the mets measured expression
mets_me <- fread(file="Z:/data/METS_model/hg19/METS_peer10_all_chr_prot_coding_gex.txt", header=T)

for (i in 1:length(mets_me$PROBE_ID)){
  mets_me$PROBE_ID[i] <- gsub('\\.[0-9]+','',mets_me$PROBE_ID[i])
}


#read in the mets predicted expression
#AFA 2 METS
afa2mets_pred <- fread(file="Z:/data/mesa_models/mets_dosages/AFA_2_METS_predicted_expression.txt", header=T)

#remove the decimal in gene id from afa2mets
for (j in 3:length(afa2mets_pred)){
  colnames(afa2mets_pred)[j] <- as.character(gsub('\\.[0-9]+','',colnames(afa2mets_pred)[j]))
}

#take the gene_ids of the genes we could predict # the gene ID no longer have extra decimal places
predcolname <- colnames(afa2mets_pred)[3:length(afa2mets_pred)]

mets_me$PROBE_ID <- as.character(mets_me$PROBE_ID)

#from the measured expression, keep only genes we could predict
mgenes_df <- subset(mets_me, mets_me$PROBE_ID %in% predcolname)

#arrange the afa2mets columns to be same order with the probe_id
afa2mets_order <- afa2mets_pred # just to make a copy of afa2mets
afa2mets_order$FID <- NULL
afa2mets_order$IID <- NULL
library(dplyr)
#selects only columns(with gene_id as colname) that is the probe_id, and keeps the new dataframe in the order of the probe_ID
afa2mets_order <- afa2mets_order %>% dplyr::select(mgenes_df$PROBE_ID)

#do the spearman correlation between the measured and predicted

afa2mets_cor <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(mgenes_df$PROBE_ID)){
  afa2mets_cor[i,1] <- mgenes_df$PROBE_ID[i]
  afa2mets_cor[i,2] <- cor(as.numeric(mgenes_df[i,2:length(mgenes_df)]), as.numeric(afa2mets_order[[i]]), method = "spearman")#calculate spearman correlation between measured and predicted gene expression
}
colnames(afa2mets_cor) <- c("gene", "spearman")
library(tidyverse)
afa2mets_cor <- drop_na(afa2mets_cor) #drop NAs
#write out results
write.table(afa2mets_cor, file = "Z:/data/mesa_models/mets_dosages/new_predixcan_AFA2METS_cor", row.names=F, quote=F, sep="\t") 



#HIS 2 METS
his2mets_pred <- fread(file="Z:/data/mesa_models/mets_dosages/HIS_2_METS_predicted_expression.txt", header=T)
#remove the decimal in gene id from afa2mets
for (j in 3:length(his2mets_pred)){
  colnames(his2mets_pred)[j] <- as.character(gsub('\\.[0-9]+','',colnames(his2mets_pred)[j]))
}

#take the gene_ids of the genes we could predict # the gene ID no longer have extra decimal places
predcolname <- colnames(his2mets_pred)[3:length(his2mets_pred)]

mets_me$PROBE_ID <- as.character(mets_me$PROBE_ID)

#from the measured expression, keep only genes we could predict
mgenes_df <- subset(mets_me, mets_me$PROBE_ID %in% predcolname)

#arrange the afa2mets columns to be same order with the probe_id
his2mets_pred$FID <-NULL
his2mets_pred$IID <- NULL
library(dplyr)
#selects only columns(with gene_id as colname) that is the probe_id, and keeps the new dataframe in the order of the probe_ID
his2mets_pred <- his2mets_pred %>% dplyr::select(mgenes_df$PROBE_ID)

#do the spearman correlation between the measured and predicted

his2mets_cor <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(mgenes_df$PROBE_ID)){
  his2mets_cor[i,1] <- mgenes_df$PROBE_ID[i]
  his2mets_cor[i,2] <- cor(as.numeric(mgenes_df[i,2:length(mgenes_df)]), as.numeric(his2mets_pred[[i]]), method = "spearman")#calculate spearman correlation between measured and predicted gene expression
}
colnames(his2mets_cor) <- c("gene", "spearman")
library(tidyverse)
his2mets_cor <- drop_na(his2mets_cor) #drop NAs
#write out results
write.table(his2mets_cor, file = "Z:/data/mesa_models/mets_dosages/new_predixcan_HIS2METS_cor", row.names=F, quote=F, sep="\t") 


#AFHI 2 METS
afhi2mets_pred <- fread(file="Z:/data/mesa_models/mets_dosages/AFHI_2_METS_predicted_expression.txt", header=T)
#remove the decimal in gene id from afa2mets
for (j in 3:length(afhi2mets_pred)){
  colnames(afhi2mets_pred)[j] <- as.character(gsub('\\.[0-9]+','',colnames(afhi2mets_pred)[j]))
}

#take the gene_ids of the genes we could predict # the gene ID no longer have extra decimal places
predcolname <- colnames(afhi2mets_pred)[3:length(afhi2mets_pred)]

mets_me$PROBE_ID <- as.character(mets_me$PROBE_ID)

#from the measured expression, keep only genes we could predict
mgenes_df <- subset(mets_me, mets_me$PROBE_ID %in% predcolname)

#arrange the afa2mets columns to be same order with the probe_id
afhi2mets_pred$FID <-NULL
afhi2mets_pred$IID <- NULL
library(dplyr)
#selects only columns(with gene_id as colname) that is the probe_id, and keeps the new dataframe in the order of the probe_ID
afhi2mets_pred <- afhi2mets_pred %>% dplyr::select(mgenes_df$PROBE_ID)

#do the spearman correlation between the measured and predicted

afhi2mets_cor <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(mgenes_df$PROBE_ID)){
  afhi2mets_cor[i,1] <- mgenes_df$PROBE_ID[i]
  afhi2mets_cor[i,2] <- cor(as.numeric(mgenes_df[i,2:length(mgenes_df)]), as.numeric(afhi2mets_pred[[i]]), method = "spearman")#calculate spearman correlation between measured and predicted gene expression
}
colnames(afhi2mets_cor) <- c("gene", "spearman")
library(tidyverse)
afhi2mets_cor <- drop_na(afhi2mets_cor) #drop NAs
#write out results
write.table(afhi2mets_cor, file = "Z:/data/mesa_models/mets_dosages/new_predixcan_AFHI2METS_cor", row.names=F, quote=F, sep="\t") 


#CAU 2 METS
cau2mets_pred <- fread(file="Z:/data/mesa_models/mets_dosages/CAU_2_METS_predicted_expression.txt", header=T)
#remove the decimal in gene id from afa2mets
for (j in 3:length(cau2mets_pred)){
  colnames(cau2mets_pred)[j] <- as.character(gsub('\\.[0-9]+','',colnames(cau2mets_pred)[j]))
}

#take the gene_ids of the genes we could predict # the gene ID no longer have extra decimal places
predcolname <- colnames(cau2mets_pred)[3:length(cau2mets_pred)]

mets_me$PROBE_ID <- as.character(mets_me$PROBE_ID)

#from the measured expression, keep only genes we could predict
mgenes_df <- subset(mets_me, mets_me$PROBE_ID %in% predcolname)

#arrange the afa2mets columns to be same order with the probe_id
cau2mets_pred$FID <-NULL
cau2mets_pred$IID <- NULL
library(dplyr)
#selects only columns(with gene_id as colname) that is the probe_id, and keeps the new dataframe in the order of the probe_ID
cau2mets_pred <- cau2mets_pred %>% dplyr::select(mgenes_df$PROBE_ID)

#do the spearman correlation between the measured and predicted

cau2mets_cor <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(mgenes_df$PROBE_ID)){
  cau2mets_cor[i,1] <- mgenes_df$PROBE_ID[i]
  cau2mets_cor[i,2] <- cor(as.numeric(mgenes_df[i,2:length(mgenes_df)]), as.numeric(cau2mets_pred[[i]]), method = "spearman")#calculate spearman correlation between measured and predicted gene expression
}
colnames(cau2mets_cor) <- c("gene", "spearman")
library(tidyverse)
cau2mets_cor <- drop_na(cau2mets_cor) #drop NAs
#write out results
write.table(cau2mets_cor, file = "Z:/data/mesa_models/mets_dosages/new_predixcan_CAU2METS_cor", row.names=F, quote=F, sep="\t")


#ALL 2 METS
all2mets_pred <- fread(file="Z:/data/mesa_models/mets_dosages/ALL_2_METS_predicted_expression.txt", header=T)
#remove the decimal in gene id from afa2mets
for (j in 3:length(all2mets_pred)){
  colnames(all2mets_pred)[j] <- as.character(gsub('\\.[0-9]+','',colnames(all2mets_pred)[j]))
}

#take the gene_ids of the genes we could predict # the gene ID no longer have extra decimal places
predcolname <- colnames(all2mets_pred)[3:length(all2mets_pred)]

mets_me$PROBE_ID <- as.character(mets_me$PROBE_ID)

#from the measured expression, keep only genes we could predict
mgenes_df <- subset(mets_me, mets_me$PROBE_ID %in% predcolname)

#arrange the afa2mets columns to be same order with the probe_id
all2mets_pred$FID <-NULL
all2mets_pred$IID <- NULL
library(dplyr)
#selects only columns(with gene_id as colname) that is the probe_id, and keeps the new dataframe in the order of the probe_ID
all2mets_pred <- all2mets_pred %>% dplyr::select(mgenes_df$PROBE_ID)

#do the spearman correlation between the measured and predicted

all2mets_cor <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(mgenes_df$PROBE_ID)){
  all2mets_cor[i,1] <- mgenes_df$PROBE_ID[i]
  all2mets_cor[i,2] <- cor(as.numeric(mgenes_df[i,2:length(mgenes_df)]), as.numeric(all2mets_pred[[i]]), method = "spearman")#calculate spearman correlation between measured and predicted gene expression
}
colnames(all2mets_cor) <- c("gene", "spearman")
library(tidyverse)
all2mets_cor <- drop_na(all2mets_cor) #drop NAs
#write out results
write.table(all2mets_cor, file = "Z:/data/mesa_models/mets_dosages/new_predixcan_ALL2METS_cor", row.names=F, quote=F, sep="\t")

