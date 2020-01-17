#make predixcan dosages with rsid

#pass arguments to the script
args <- commandArgs(trailingOnly = T)
chrom <- as.character(args[1])

library(data.table)
library(dplyr)
print(chrom)

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

"%&%" <- function(a,b) paste(a,b, sep = "")

caupheno <- read.table(file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/chr" %&% chrom %&% "txt.gz", header=T)
for (j in 7:length(caupheno)){
  colnames(caupheno)[j] <- strsplit(colnames(caupheno)[j], "_")[[1]][2]
}# just to convert the columns to only string of numbers

caupheno <- caupheno[, (names(caupheno) %in% header)] #keep only columns that are in the header
#Now write the dosage file back
caupheno$snp_ID <- as.character(caupheno$snp_ID)

#change the dosages to shape for doing predixcan hg19 b37
for (k in 1:nrow(caupheno)){
  caupheno$snp_ID[k] <- as.character(caupheno$chr[k] %&% "_" %&% caupheno$pos[k] %&% "_" %&% caupheno$ref[k] %&% "_" %&% caupheno$alt[k] %&% "_b37")
  caupheno$snp_ID <- as.character(caupheno$snp_ID)
}

#check the weight file, and use it to change to dosages to have same rsid
allwei <- fread(file="/home/pokoro/MESA.ALL.WG.PC3.PF10.unpruned.rsid.hg19.weights.txt", header=F)
allwei$V3 <- as.character(allwei$V3)
allwei[1:5,1:6]

rsid2 <- inner_join(allwei, caupheno, by = c("V3"="snp_ID"))
rsid2 <- rsid2[,c(2,7:length(rsid2))] #take out rsid and the other dosage columns ie chr pos ref alt maf id1....idn
rsid2 <- rsid2[,c(2,1,3:length(rsid2))]
rsid2 <- unique(rsid2) #remove duplicates
rsid2[1:5,1:10]


#write out without column names. chromosome rsid position allele1 allele2 MAF id1 ..... idn. predixcan does not need colnames
write.table(rsid2, file = "/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/predixcan/correct_dosages/chr" %&% chrom %&% "_cau.txt", quote=F, row.names=F, sep ="\t", col.names=F)


