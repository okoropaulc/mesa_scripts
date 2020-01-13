#association test of ml predicted expressions vs phenotypes
#the predicted expression file format should be like this
#sample_id gene1 gene2 gene3 #that is the column


library(data.table)

#cauphenoall <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr1_all_sample.txt", header=T, nrows=3)
#cau_svr <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_8.txt", header=T)
#cau_svr14 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_14.txt", header=T)
#cau_pd <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/grid_optimized_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_predicted_gene_expr_chr8.txt",header=T)
#cau_pd14 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/grid_optimized_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_predicted_gene_expr_chr14.txt",header=T)

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

#change the sample ids to numbers only
#for (i in 1:nrow(cau_pd)){
#  cau_pd$V1[i] <- strsplit(cau_pd$V1[i], "X")[[1]][2]
#}


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

# gene <- genes[1]
# pred_gene_exp <- rfmerged[[gene]]
# model <- lm(phenotype ~ pred_gene_exp, data = rfmerged)
# results <- coef(summary(model))[c(2,6,8,4)]

#line <- c(gene,results)
#assoc_df <- NULL
#assoc_df <- rbind(assoc_df,line)
#colnames(assoc_df) <- c("gene", "beta", "t", "p", "se(beta)")

#assoc_df <-as.data.frame(assoc_df)

#merged2 <- merged1

#merged2$sidno <- NULL
#colnames(merged2)[1] <- "phenotype"#rank_plt5
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

#check gene name
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
#fwrite(rfassoc, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/rf_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5_1207.txt",row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/svr_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5_1207.txt",row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/knn_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5_1207.txt",row.names=F, quote=F, sep="\t")


#write_association <- function(assoc_df, output_file) {
#  write.table(assoc_df, output_file, col.names = T, row.names = F, quote = F)
#}

rfassoc1 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/rf_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt", header=T)
svrassoc1 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/svr_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt", header=T)
knnassoc1 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/knn_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5.txt", header=T)

rfassoc <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/rf_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5_1207.txt", header=T)
svrassoc <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/svr_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5_1207.txt", header=T)
knnassoc <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/knn_pred_ex_pheno_assoc_all_2_cau_thrombo_rankplt5_1207.txt", header=T)


#make qqplot
#install.packages("qqman")
library(qqman)

qq(knnassoc1$p, main="knn")
qq(svrassoc1$p, main="svr")
qq(rfassoc1$p, main="rf")

#qqnorm(knnassoc1$p) shows the pvalue came from a population that is normal distribution


library(data.table)
library(dplyr)
#filter the pheno asscoiation result to have only ALL gene models with cv R2 > 0.01
#read in the best gird results for the ALL
all_rfgrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_rf_all_chrom.txt", header=T)
all_svrgrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_svr_all_chrom.txt", header=T)
all_knngrid <- fread(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/ALL_best_grid_knn_all_chrom.txt", header=T)

all_rfgrid$Gene_ID <- as.character(all_rfgrid$Gene_ID)
all_svrgrid$Gene_ID <- as.character(all_svrgrid$Gene_ID)
all_knngrid$Gene_ID <- as.character(all_knngrid$Gene_ID)

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

rfassoc_cvr0.01 <- inner_join(all_rfgrid, rfassoc, by = c("Gene_ID" = "gene"))
svrassoc_cvr0.01 <- inner_join(all_svrgrid, svrassoc, by = c("Gene_ID" = "gene"))
knnassoc_cvr0.01 <- inner_join(all_knngrid, knnassoc, by = c("Gene_ID" = "gene"))

#qq(rfassoc_cvr0.01$p, main="rf0.01")
