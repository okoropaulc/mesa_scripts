#pass arguments to the script
args <- commandArgs(trailingOnly = T)
chrom <- as.character(args[1])

library(data.table)
print(chrom)

#mycau <- read.table(file="Z:/data/mesa_models/mesa_pheno/thrombotic/predixcan/chr22_cau.txt", nrows=3)
#mycau2 <- read.table(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/cau_dosage_chr1_all_sample.txt", header=3, nrows=3)

#pcau <- read.table(file="Z:/data/chr22.txt.gz", nrows=3) #pfiorica dosage

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


"%&%" <- function(a,b) paste(a,b, sep = "")

caupheno <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/chr" %&% chrom %&% "txt.gz", header=T)
for (j in 7:length(caupheno)){
  colnames(caupheno)[j] <- strsplit(colnames(caupheno)[j], "_")[[1]][2]
}# just to convert the columns to only string of numbers
  
caupheno <- caupheno[, (names(caupheno) %in% header)] #keep only columns that are in the header
#Now write the dosage file back
caupheno$snp_ID <- as.character(caupheno$snp_ID)
  
#change the dosages to shape for doing predixcan hg38
for (k in 1:nrow(caupheno)){
  caupheno$snp_ID[k] <- as.character("hg38:" %&% caupheno$chr[k] %&% ":" %&% caupheno$pos[k])
  caupheno$snp_ID <- as.character(caupheno$snp_ID)
}

#write out without column names. chromosome rsid position allele1 allele2 MAF id1 ..... idn. predixcan does not need colnames
write.table(caupheno, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan/chr" %&% chrom %&% "_cau.txt", quote=F, row.names=F, sep ="\t", col.names=F)


# for (i in 1:22){
#   print(i)
#   
#   caupheno <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/chr" %&% i %&% "txt.gz", header=T)
#   for (j in 7:length(caupheno)){
#     colnames(caupheno)[j] <- strsplit(colnames(caupheno)[j], "_")[[1]][2]
#   }# just to convert the columns to only string of numbers
#   
#   caupheno <- caupheno[, (names(caupheno) %in% header)] #keep only columns that are in the header
#   #Now write the dosage file back
#   caupheno$snp_ID <- as.character(caupheno$snp_ID)
#   
#   #change the dosages to shape for doing predixcan hg38
#   for (k in 1:nrow(caupheno)){
#     caupheno$snp_ID[k] <- as.character("hg38:" %&% caupheno$chr[k] %&% ":" %&% caupheno$pos[k])
#     caupheno$snp_ID <- as.character(caupheno$snp_ID)
#   }
#   #write out without column names. chromosome rsid position allele1 allele2 MAF id1 ..... idn
#   write.table(caupheno, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan/chr" %&% i %&% "_cau.txt", quote=F, row.names=F, sep ="\t", col.names=F)
#   
# }