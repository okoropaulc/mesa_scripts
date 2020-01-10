#lipid_pheno <- read.csv(file = "Z:/data/mesa_models/mesa_pheno/Exam1Main.csv", header = T)

thrombomodulin <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/MESA_thrombotic_phenotypes.txt", header=T)

#check CAU dosages for thrombomodulin
print("check CAU dosages for thrombomodulin")

caudos1 <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/chr1txt.gz", header=T, nrows = 3)
header3 <- colnames(caudos1)[1:6] #store the names of the first six column

caudos1[,c(1:6)] <- NULL



#strsplit returns list, and to access it, use [[]] to access a row and [] to access a column

#change the colnames of cau dosages to be numbers only
#print("change the colnames of cau dosages to be numbers only")

for (i in 1:length(caudos1)){
  colnames(caudos1)[i] <- strsplit(colnames(caudos1)[i], "_")[[1]][2]
}

#read in the mesa cau samples used to build the model
mcaudos1 <- read.table(file = "/home/pokoro/data/mesa_models/cau/CAU_1_snp.txt", header=T, nrow=1)
mcaudos1$id <- NULL

#change the colnames to numbers only
#print("change the colnames to numbers only")
for (i in 1:length(mcaudos1)){
  colnames(mcaudos1)[i] <- strsplit(colnames(mcaudos1)[i], "X")[[1]][2]
}

#read in the columns
#print("read in the columns")
caucol <- colnames(caudos1)
caucol <- sort(caucol)

mcaucol <- colnames(mcaudos1)
mcaucol <- sort(mcaucol)


#print("keep sample ID's that are not in the dosages used to build model")
fcau <- caudos1[, !(names(caudos1) %in% mcaucol)] #keep sample ID's that are not in the dosages used to build model


#drop NA in phenotype
#print("drop NA in phenotype")
library(tidyverse)

#print("drop NA")
#keep ID and rankplt5
thrombomodulin <- thrombomodulin[,c(1,5)]
thrombomodulin <- drop_na(thrombomodulin) #drop NA

#Now keep only sample ID in dosages that are in phenotype
#print("Now keep only sample ID in dosages that are in phenotype")
fcau <- fcau[, (names(fcau) %in% thrombomodulin$sidno)]

print("rename the row with the sample ID")
rownames(thrombomodulin) <- thrombomodulin$sidno #rename the row with the sample ID

#Drop rows where sample ID is not the phenotype dosages

remthromb <- thrombomodulin[(rownames(thrombomodulin) %in% names(fcau)),]

#Therefore, number of samples we have thrombomoodulin phenotype and dosages not in training data is 665

remid <- names(fcau) #write the remaining sample IDs

header <- c(header3, remid) #combine the first 6 column with the remaining sample ID

#Now read in the genotype files and select only the sample IDs with phenotype

library(data.table)

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
  
  write.table(caupheno, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr" %&% i %&% ".txt", quote=F, row.names=F, sep ="\t")
  
}



#change the dosages to shape for doing gene expression imputation for build 37
#chr_pos_ref_alt_build

"%&%" <- function(a,b) paste(a,b, sep = "")

for (j in 1:22){
  
  caupheno <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr" %&% j %&% ".txt", header=T)
  
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
  write.table(cau_dosage, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_nocautrain_chr" %&% j %&% ".txt", quote=F, row.names=F, sep ="\t")
  
}


#rfgex <- read.table(file="Z:/data/mesa_models/python_ml_models/results/CAU_2_METS_rf_predicted_gene_expr_chr1.txt", header=T)

#caupheno <- read.table(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr1_all_sample.txt", header=T, nrow=3)
#mcaudos1 <- read.table(file = "Z:/data/mesa_models/cau/CAU_1_snp.txt", header=T, nrow=3)
#caudoschr <- read.table(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/chr1txt.gz", header=T, nrows = 3)

#for (i in 6:length(caudoschr)){
#  colnames(caudoschr)[i] <- strsplit(colnames(caudoschr)[i], "_")[[1]][2]
#}

#caudoschr <- caudoschr[, (names(caudoschr) %in% header)]

#mcaudos1 <- read.table(file = "Z:/data/mesa_models/cau/CAU_1_snp.txt", header=T, nrow=1)

#write.table(caudoschr, file = "Z:/data/mesa_models/mesa_pheno/thrombotic/try2.txt", quote=F, sep="\t", row.names=F)

#sidall <- read.table(file="Z:/data/mesa_models/all/whole_genotypes/ALL.chr1.genotype.txt.gz", header=T, nrows=3)
#mcaudos1 <- read.table(file = "Z:/data/mesa_models/cau/CAU_1_snp.txt", header=T, nrow=3)
#mhisdos1 <- read.table(file = "Z:/data/mesa_models/his/HIS_1_snp.txt", header=T, nrow=3)
#mafados1 <- read.table(file = "Z:/data/mesa_models/AFA_1_snp.txt", header=T, nrow=3)
