#Make snp annotation file for the cau pheno dosages

# read in example snp annotation of cau
#mcauannot <- read.table(file = "Z:/data/mesa_models/cau/CAU_1_annot.txt", header=T, nrow=3)

#pheno cau dosage
#caudos1 <- read.table(file = "Z:/data/mesa_models/mesa_pheno/thrombotic/chr1txt.gz", header=T, nrows = 3)


"%&%" <- function(a,b) paste(a,b, sep = "")

library(data.table)
#build snp annot file for each of the chrom


print("loop started")

for (j in 1:22){
  print(j)
  #cau pheno imputation dosage
  cau_imp <- fread("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr" %&% j %&%".txt", header=T, select="id", showProgress=T)
  cau_imp$id <- as.character(cau_imp$id)
  
  
  #create the annot df
  snpannot <- data.frame(chr="",pos="",varID="",refAllele="",effectAllele="",rsid="", stringsAsFactors=F)
  
  for (i in 1:nrow(cau_imp)){
    #break the id by "_"
    id_split <- strsplit(cau_imp$id[i], "_")
    snpannot[i,] <- c(id_split[[1]][1], id_split[[1]][2], cau_imp$id[i], id_split[[1]][3], id_split[[1]][4], 
                      id_split[[1]][1] %&% ":" %&% id_split[[1]][2])
  }
  fwrite(snpannot, file="/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr" %&% j %&% "_annot.txt", 
              row.names=F, quote=F, sep="\t", showProgress = T)
  print("done")
}
