
#strsplit returns list, and to access it, use [[]] to access a row and [] to access a column

#change the dosages to shap for doing gene expression imputation for build 37
#chr_pos_ref_alt_build

#pass arguments to the script
args <- commandArgs(trailingOnly = T)
j <- as.character(args[1])
pop <- as.character(args[2])

"%&%" <- function(a,b) paste(a,b, sep = "")

dosage <- read.table(file = "/home/pokoro/data/lauren_mesa/"%&%pop%&%"_dosages/chr" %&% j %&% ".maf0.01.r20.8noambig.dosage.txt.gz", header=F)

dos_id <- data.frame(id=rep("", nrow(dosage))) #create an "id" column with NA's

join_dosage <- cbind(dos_id, dosage) #join the id with the dosage dataframe

#make the id to be in this format #chr_pos_ref_alt_build

join_dosage$id <- as.character(join_dosage$id)
#dosage$chr <- as.character(dosage$chr) #chromosome number. I commented it out because I just want to use the input chrom no
dosage$V3 <- as.character(dosage$V3) #snp position pos
dosage$V4 <- as.character(dosage$V4) #ref allele
dosage$V5 <- as.character(dosage$V5) #alt allele

for (i in 1:nrow(dosage)){
  join_dosage$id[i] <- as.character(j) %&% "_" %&% dosage$V3[i] %&% "_" %&% dosage$V4[i] %&% "_" %&% dosage$V5[i] %&% "_b37"
  join_dosage$id <- as.character(join_dosage$id)
}

join_dosage[,c(2:7)] <- NULL #remove the unwanted columns, that is chr, rsid, pos, ref, alt, aa_freq
join_dosage$id <- as.character(join_dosage$id)
#write out the file
write.table(join_dosage, file = "/home/pokoro/data/lauren_mesa/ml_dosages/"%&%pop%&%"/chr" %&% j %&% ".txt", quote=F, row.names=F, sep ="\t")
