#lipid_pheno <- read.csv(file = "Z:/data/mesa_models/mesa_pheno/Exam1Main.csv", header = T)

thrombomodulin <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)

#check CAU dosages for thrombomodulin
#print("check CAU dosages for thrombomodulin")

caudos1 <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/chr1txt.gz", header=T, nrows = 3)
header3 <- colnames(caudos1)[1:6] #store the names of the first six column

caudos1[,c(1:6)] <- NULL



#strsplit returns list, and to access it, use [[]] to access a row and [] to access a column

#change the colnames of cau dosages to be numbers only
#print("change the colnames of cau dosages to be numbers only")

for (i in 1:length(caudos1)){
  colnames(caudos1)[i] <- strsplit(colnames(caudos1)[i], "_")[[1]][2]
}

#read in the columns
#print("read in the columns")
caucol <- colnames(caudos1)
caucol <- sort(caucol)


#drop NA in phenotype
#print("drop NA in phenotype")
library(tidyverse)

#print("drop NA")
#keep ID and rankplt5
thrombomodulin <- thrombomodulin[,c(1,5)]
thrombomodulin <- drop_na(thrombomodulin) #drop NA

#Now keep only sample ID in dosages that are in phenotype
#print("Now keep only sample ID in dosages that are in phenotype") 2750 samples
fcau <- thrombomodulin[,1]

#print("rename the row with the sample ID")
rownames(thrombomodulin) <- thrombomodulin$sidno #rename the row with the sample ID


header <- c(header3, fcau) #combine the first 6 column with the remaining sample ID

#Now read in the genotype files and select only the sample IDs with phenotype

#library(data.table)

#use loop to go through all the 22 chromosomes
"%&%" <- function(a,b) paste(a,b, sep = "")

for (i in 1:22){
  print(i)
  
  caupheno <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/chr" %&% i %&% "txt.gz", header=T)
  for (j in 6:length(caupheno)){
    colnames(caupheno)[j] <- strsplit(colnames(caupheno)[j], "_")[[1]][2]
  }# just to convert the columns to only string of numbers
  
  caupheno <- caupheno[, (names(caupheno) %in% header)] #keep only columns that are in the header
  #Now write the dosage file back
  
  write.table(caupheno, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr" %&% i %&% "_all_sample.txt", quote=F, row.names=F, sep ="\t")
  
}


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
