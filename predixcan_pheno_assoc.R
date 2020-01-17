

library(data.table)
#this is the old one with zrosfread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/en_cau_rankplt5_predicted_expression.txt", header=T)
pxcan <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/en_cau_rankplt5_predicted_expression.txt", header=T)
#header <- colnames(mycau)[1:5]
#mycau[1:7,1:5]
#predrf <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/pred_expr/chunk/transformed_full_chrom_ALL_2_CAU_thrombomodulin_rankplt5_pheno_rf_pred_expr.txt",header=T)

thrombomodulin <- fread(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)
library(tidyverse)

#drop NA
#keep ID and rankplt5
thrombomodulin <- thrombomodulin[,c(1,5)]
thrombomodulin <- drop_na(thrombomodulin) #drop NA
thrombomodulin$sidno <- as.character(thrombomodulin$sidno)
colnames(thrombomodulin)[2] <- "phenotype" #rank_plt5

#pxcan1 <- drop_na(pxcan)
#remove one of the ID
pxcan$FID <- NULL


pxcan$IID <- as.character(pxcan$IID)

library(dplyr)

#merge pheno and pred expr on their common sample
#take the genes
pxcangenes <- colnames(pxcan)[2:length(pxcan)]

pxcanmerged <- inner_join(thrombomodulin, pxcan, by = c("sidno" = "IID"))

pxcanmerged$sidno <- NULL

#functionize #pheno assoc function
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

pxcanassoc <- association(pxcanmerged, pxcangenes)
#write the assoc result
fwrite(pxcanassoc, file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/en_cau_rankplt5_association.txt", quote=F, row.names=F, sep ="\t")

#remove NAs
pxcanassoc1 <- drop_na(pxcanassoc)

#check gene name
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)

#check why we get NA for some gene
#gene TRPM4 ENSG00000130529.11
pred_gene_exp <- pxcanmerged[["ENSG00000130529.11"]]
model <- lm(phenotype ~ pred_gene_exp, data = pxcanmerged)
results <- coef(summary(model))[c(2,6,8,4)] #therefore the reason for the NA is because the predicted expression for the gene were 0
#predicted expression for all genes were zero

#the reason could be because of no overlap between the snps in the test cohort dosages, and the ALL model
#so lets check
#first, read in the file ryan used to do the liftover from hg19 to cpos
cpos <- fread(file="Z:/ALL_hg38_cpos_varid.txt", header=F)
cpos[1:3,] #just to see the file header, cause opening the whole file crashes rstudio
# 3rd colum is cpos
cpos$V3 <- as.character(cpos$V3)
allcpos <- as.data.frame(cpos$V3)
colnames(allcpos) <- "cpos"
allcpos$cpos <- as.character(allcpos$cpos)
#read in chrom 22 dosages
chr22 <- fread(file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/chr21_cau.txt.gz", header=F)
chr22[1:3,1:3]#just to see the header
#2nd column is the cpos
chr22$V2 <- as.character(chr22$V2)
chr22cpos <- as.data.frame(chr22$V2)
colnames(chr22cpos) <- "cpos"
chr22cpos$cpos <- as.character(chr22cpos$cpos)
#now find their overlap
library(dplyr)
chroverlap <- inner_join(allcpos, chr22cpos, by = c("cpos"="cpos"))
#In conclusion, the dosage file is bad