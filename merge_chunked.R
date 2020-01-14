#merge all the chunked expression prediction results
"%&%" <- function(a,b) paste(a,b, sep = "")
library(data.table)

algs <- c("rf", "svr", "knn")
rf <- NULL
svr <- NULL
knn <- NULL

#eg file; chr9_chunk5_ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr.txt
#chunked chroms are
chunked <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,19,20,22)
for (i in chunked){
  rf1 <- NULL
  svr1 <- NULL
  knn1 <- NULL
  for (j in 1:5){
    rf1 <- rbind(rf1, fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr" %&% i %&% "_chunk" %&% j %&% "_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt", header=T))
    svr1 <- rbind(svr1, fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr" %&% i %&% "_chunk" %&% j %&% "_ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr.txt", header=T))
    knn1 <- rbind(knn1, fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr" %&% i %&% "_chunk" %&% j %&% "_ALL_2_CAU_thrombomodulin_rankplt5_phenoknn_pred_expr.txt", header=T))
    rf <- rbind(rf, fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr" %&% i %&% "_chunk" %&% j %&% "_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt", header=T))
    svr <- rbind(svr, fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr" %&% i %&% "_chunk" %&% j %&% "_ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr.txt", header=T))
    knn <- rbind(knn, fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr" %&% i %&% "_chunk" %&% j %&% "_ALL_2_CAU_thrombomodulin_rankplt5_phenoknn_pred_expr.txt", header=T))
    
  }
  fwrite(rf1, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/full_chr" %&% i %&% "_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_expr.txt", row.names=F, quote=F, sep="\t")
  fwrite(svr1, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/full_chr" %&% i %&% "_ALL_2_CAU_thrombomodulin_rankplt5_pheno_svr_pred_expr.txt", row.names=F, quote=F, sep="\t")
  fwrite(knn1, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/full_chr" %&% i %&% "_ALL_2_CAU_thrombomodulin_rankplt5_pheno_knn_pred_expr.txt", row.names=F, quote=F, sep="\t")
}


#read in the other chroms that were completed before the genotype was chunked chroms 13,18, 21
chr13knn <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenoknn_pred_expr_13.txt", header=T)
chr13svr <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_13.txt", header=T)
chr13rf <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr_13.txt", header=T)

chr18knn <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenoknn_pred_expr_18.txt", header=T)
chr18svr <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_18.txt", header=T)
chr18rf <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr_18.txt", header=T)

chr21knn <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenoknn_pred_expr_21.txt", header=T)
chr21svr <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_21.txt", header=T)
chr21rf <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr_21.txt", header=T)
#merge them with the all
rf <- rbind(rf, chr21rf, chr18rf, chr13rf)
svr <- rbind(svr, chr21svr, chr18svr, chr13svr)
knn <- rbind(knn, chr21knn, chr18knn, chr13knn)


fwrite(rf, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_expr.txt", row.names=F, quote=F, sep="\t")
fwrite(svr, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_svr_pred_expr.txt", row.names=F, quote=F, sep="\t")
fwrite(knn, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_knn_pred_expr.txt", row.names=F, quote=F, sep="\t")


#transform the pred exp to the required
#rf
trf <- as.data.frame(t(rf))
rowrf <- rf$gene_id #take the gene id's
colnames(trf) <- rowrf
sampleidrf <- colnames(rf)[2:length(rf)] #take sample only ids, the number is constant
#change the sample ids to numbers only
for (i in 1:length(sampleidrf)){
  sampleidrf[i] <- strsplit(sampleidrf[i], "X")[[1]][2]
}
trf <- trf[-1,] #delete the first row
trf <- cbind(sampleidrf, trf)

#svr
tsvr <- as.data.frame(t(svr))
rowsvr <- svr$gene_id #take the gene id's
colnames(tsvr) <- rowsvr
sampleidsvr <- colnames(svr)[2:length(svr)] #take sample only ids, the number is constant
#change the sample ids to numbers only
for (i in 1:length(sampleidsvr)){
  sampleidsvr[i] <- strsplit(sampleidsvr[i], "X")[[1]][2]
}
tsvr <- tsvr[-1,] #delete the first row
tsvr <- cbind(sampleidsvr, tsvr)

#knn
tknn <- as.data.frame(t(knn))
rowknn <- knn$gene_id #take the gene id's
colnames(tknn) <- rowknn
sampleidknn <- colnames(knn)[2:length(knn)] #take sample only ids, the number is constant
#change the sample ids to numbers only
for (i in 1:length(sampleidknn)){
  sampleidknn[i] <- strsplit(sampleidknn[i], "X")[[1]][2]
}
tknn <- tknn[-1,] #delete the first row
tknn <- cbind(sampleidknn, tknn)

fwrite(trf, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_expr.txt", row.names=F, quote=F, sep="\t")
fwrite(tsvr, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_svr_pred_expr.txt", row.names=F, quote=F, sep="\t")
fwrite(tknn, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_knn_pred_expr.txt", row.names=F, quote=F, sep="\t")
