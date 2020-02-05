#pass arguments to the script
args <- commandArgs(trailingOnly = T)
j <- as.character(args[1])
pop <- as.character(args[2])

"%&%" <- function(a,b) paste(a,b, sep = "")

library(data.table)
#build snp annot file for each of the chrom
print(j)

#cau pheno imputation dosage
mesa <- fread("/home/pokoro/data/lauren_mesa/ml_dosages/"%&%pop%&%"/chr" %&% j %&%".txt", header=T, select="id")
mesa$id <- as.character(mesa$id)


#create the annot df
snpannot <- data.frame(chr="",pos="",varID="",refAllele="",effectAllele="",rsid="", stringsAsFactors=F)

for (i in 1:nrow(mesa)){
  #break the id by "_"
  id_split <- strsplit(mesa$id[i], "_")
  #the snpannot[1,] is to store only one row in the df, and save memory
  snpannot[1,] <- c(id_split[[1]][1], id_split[[1]][2], mesa$id[i], id_split[[1]][3], id_split[[1]][4], 
                    id_split[[1]][1] %&% ":" %&% id_split[[1]][2])
  fwrite(snpannot, file="/home/pokoro/data/lauren_mesa/snp_annotation/"%&%pop%&%"/chr" %&% j %&% "_annot.txt", 
         row.names=F, quote=F, sep="\t", append = T)
}
