#change the dosages to shap for doing gene expression imputation for build 37
#chr_pos_ref_alt_build

"%&%" <- function(a,b) paste(a,b, sep = "")

for (j in 1:22){
  
  caupheno <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr" %&% j %&% "_all_sample.txt", header=T)
  
  cau_id <- data.frame(id=rep("", nrow(caupheno))) #create an "id" column with NA's
  
  cau_dosage <- cbind(cau_id, caupheno) #join the id with the dosage dataframe
  
  #make the id to be in this format #chr_pos_ref_alt_build
  
  cau_dosage$id <- as.character(cau_dosage$id)
  caupheno$chr <- as.character(caupheno$chr)
  caupheno$pos <- as.character(caupheno$pos)
  caupheno$ref <- as.character(caupheno$ref)
  caupheno$alt <- as.character(caupheno$alt)
  
  for (i in 1:nrow(caupheno)){
    cau_dosage$id[i] <- as.character(caupheno$chr[i]) %&% "_" %&% caupheno$pos[i] %&% "_" %&% caupheno$ref[i] %&% "_" %&% caupheno$alt[i] %&% "_b37"
    cau_dosage$id <- as.character(cau_dosage$id)
  }
  
  cau_dosage[,c(2:6)] <- NULL #remove the unwanted columns
  cau_dosage$id <- as.character(cau_dosage$id)
  #write out the file
  write.table(cau_dosage, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr" %&% j %&% ".txt", quote=F, row.names=F, sep ="\t")
  
}

#cau_imp <- read.table("Z:/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr1.txt", header=T, nrows=3)
