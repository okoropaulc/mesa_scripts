#association test of ml predicted expressions vs phenotypes
#the predicted expression file format should be like this
#sample_id gene1 gene2 gene3 #that is the column


library(data.table)

#read in pred expr
predrf <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_expr.txt",header=T)
predsvr <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_svr_pred_expr.txt",header=T)
predknn <- fread(file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_knn_pred_expr.txt",header=T)

#
thrombomodulin <- fread(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)
library(tidyverse)

#drop NA
#keep ID and rankplt5
thrombomodulin <- thrombomodulin[,c(1,5)]
thrombomodulin <- drop_na(thrombomodulin) #drop NA
thrombomodulin$sidno <- as.character(thrombomodulin$sidno)
colnames(thrombomodulin)[2] <- "phenotype" #rank_plt5

predrf$sampleidrf <- as.character(predrf$sampleidrf)
predsvr$sampleidsvr <- as.character(predsvr$sampleidsvr)
predknn$sampleidknn <- as.character(predknn$sampleidknn)


library(dplyr)

#merge pheno and pred expr on their common sample
rfgenes <- colnames(predrf)[2:length(predrf)]
svrgenes <- colnames(predsvr)[2:length(predsvr)]
knngenes <- colnames(predknn)[2:length(predknn)]


rfmerged <- inner_join(thrombomodulin, predrf, by = c("sidno" = "sampleidrf"))
svrmerged <- inner_join(thrombomodulin, predsvr, by = c("sidno" = "sampleidsvr"))
knnmerged <- inner_join(thrombomodulin, predknn, by = c("sidno" = "sampleidknn"))

rfmerged$sidno <- NULL
svrmerged$sidno <- NULL
knnmerged$sidno <- NULL

# functionize


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

rfassoc <- association(rfmerged, rfgenes)
svrassoc <- association(svrmerged, svrgenes)
knnassoc <- association(knnmerged, knngenes)

#write out the result
fwrite(rfassoc, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/rf_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt",row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/svr_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt",row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/knn_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt",row.names=F, quote=F, sep="\t")

