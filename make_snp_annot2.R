#Make snp annotation file for the cau pheno dosages

# read in example snp annotation of cau
#mcauannot <- read.table(file = "Z:/data/mesa_models/cau/CAU_1_annot.txt", header=T, nrow=3)

#pheno cau dosage
#caudos1 <- read.table(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/chr1txt.gz", header=T, nrows = 3)

#pass arguments to the script
args <- commandArgs(trailingOnly = T)
j <- as.character(args[1])


"%&%" <- function(a,b) paste(a,b, sep = "")

library(data.table)
#build snp annot file for each of the chrom
print(j)

#cau pheno imputation dosage
cau_imp <- fread("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr" %&% j %&%".txt", 
                 header=T, select="id", showProgress=T)
cau_imp$id <- as.character(cau_imp$id)
  
  
#create the annot df
snpannot <- data.frame(chr="",pos="",varID="",refAllele="",effectAllele="",rsid="", stringsAsFactors=F)
  
for (i in 1:nrow(cau_imp)){
  #break the id by "_"
  id_split <- strsplit(cau_imp$id[i], "_")
  #the snpannot[1,] is to store only one row in the df, and save memory
  snpannot[1,] <- c(id_split[[1]][1], id_split[[1]][2], cau_imp$id[i], id_split[[1]][3], id_split[[1]][4], 
                    id_split[[1]][1] %&% ":" %&% id_split[[1]][2])
  fwrite(snpannot, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr" %&% j %&% "_annot.txt", 
         row.names=F, quote=F, sep="\t", append = T)
}


#cau_desk <- fread(file="C:/Users/okoro/OneDrive/Desktop/cau_imputation_dosage_chr" %&% j %&% "_annot.txt", header=T)
#cau_desk <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr1_annot.txt", header=T)

#cau_desk <- fread(file="C:/Users/okoro/OneDrive/Desktop/gex.txt", header=T)
#cau_rf <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr_21.txt", header=T)
#cau_svr <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_22.txt", header=T)
#cau_knn <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenoknn_pred_expr_22.txt", header=T)
#cau_knn <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr8_chunk1_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_gene_expr.txt",header=T)
cau_svr <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_8.txt", header=T)
cau_svr14 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/ALL_2_CAU_thrombomodulin_rankplt5_phenosvr_pred_expr_14.txt", header=T)

cau_rf81 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr19_chunk1_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt",header=T)
cau_rf82 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr19_chunk2_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt",header=T)
cau_rf83 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr19_chunk3_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt",header=T)
cau_rf84 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr19_chunk4_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt",header=T)
cau_rf85 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr19_chunk5_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt",header=T)

cau_rf145 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/chr14_chunk5_ALL_2_CAU_thrombomodulin_rankplt5_phenorf_pred_expr.txt",header=T)

library(data.table)

cau_pd <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/grid_optimized_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_predicted_gene_expr_chr8.txt",header=T)
cau_pd14 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/grid_optimized_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_predicted_gene_expr_chr14.txt",header=T)
