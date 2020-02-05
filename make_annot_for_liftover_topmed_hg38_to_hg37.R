#Make snp annotation file for the topmed dosages for liftover purposes

#format is
#chr_no     snp_pos     snp_pos+1     snp_ID(use one in the dosage file)

#AFA

#pass arguments to the script
args <- commandArgs(trailingOnly = T)
j <- as.character(args[1])
pop <- as.character(args[2])

"%&%" <- function(a,b) paste(a,b, sep = "")

library(data.table)
#build snp annot file for each of the chrom
print(j)

#cau pheno imputation dosage
topmed <- fread("/home/pokoro/data/mesa_models/topmed_dosages/" %&% pop %&% "/chr" %&% j %&%".maf0.01.R20.8.dosage.txt.gz", 
                 header=T, select=c("chr","snp_ID","pos"), showProgress=T)



#create the annot df
snpannot <- data.frame(chr="",pos="",pos1="",snp_ID="", stringsAsFactors=F)

for (i in 1:nrow(topmed)){
  #the snpannot[1,] is to store only one row in the df, and save memory
  snpannot[1,] <- c(topmed$chr[i],topmed$pos[i],(topmed$pos[i]+1),topmed$snp_ID[i])
  fwrite(snpannot, file="/home/pokoro/data/mesa_models/topmed_dosages/" %&% pop %&% "/chr" %&% j %&% "_annot.txt", 
         row.names=F, quote=F, sep="\t", append = T)
}

#chr22 <- fread(file="Z:/data/mesa_models/topmed_dosages/AFA/chr22_annot.txt", header=T)
