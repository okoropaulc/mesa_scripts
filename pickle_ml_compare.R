library(data.table)
library(dplyr)
library(ggplot2)

"%&%" <- function(a,b) paste(a,b, sep = "")
pop <- "his"

rf <- fread(file = "Z:/data/twas_mesa/ml_pred/transformed_" %&% tolower(pop) %&% "_allchrom_r2_0.01_rf.txt", header = T)

# rf_em <- fread(file="Z:/data/ml_software/" %&% pop %&% "_RF_missing_snps_imputed_with_EM_predicted_expression_chrom_21_transposed.txt",
#                header=T)
rf_em <- NULL
#combine the pickle rf into all chrom
for (i in 1:22){
  rf_em <- cbind(rf_em, fread(file="Z:/data/ml_software/" %&% pop %&% "_RF_missing_snps_imputed_with_EM_predicted_expression_chrom_" %&% i %&% "_transposed.txt",
                              header=T))
}

#write out the whole chrom
fwrite(rf_em, file="Z:/data/ml_software/" %&% pop %&% "_RF_missing_snps_imputed_with_EM_predicted_expression_allchrom_transposed.txt",
       row.names=F, quote=F, sep="\t")


em <- rf_em[,2:length(rf_em)]

#remove the decimals in the gene id (colnames) of the predicted expression
for (i in 1:ncol(em)){
  colnames(em)[i] <- gsub('\\.[0-9]+','',colnames(em)[i])
}

#take only the chr22 genes in rf_em
#chr22rf <- rf %>% select(colnames(em))


old_rf <- rf[,2:length(rf)]

#compare the means of each gene prediction between old rf and pickle rf which involved imputing missing snps with EM

old_rf_mean <- as.data.frame(colMeans(old_rf)) #as.data.frame(colMeans(em[,sapply(em, as.numeric)]))#
old_rf_mean$gene <- colnames(old_rf)
names(old_rf_mean) <- c("meanrf","gene")

em_mean <- as.data.frame(colMeans(em[,sapply(em, as.numeric)]))#as.data.frame(colMeans(em))
em_mean$gene <- colnames(em)
names(em_mean) <- c("meanem","gene")

merged <- inner_join(old_rf_mean, em_mean, by = c("gene"="gene"))


tiff("Z:/data/ml_software/his_rf_pickle_compare.tiff", width = 7, height = 5, units = 'cm', res = 300, compression = 'lzw')

ggplot(merged, aes(x=meanrf, y=meanem)) + geom_point() + theme_classic(16) + xlab("No Missing SNPs") + 
  ylab("Missing SNPs Imputed with EM")

dev.off()

#aa <- as.data.frame(colMeans(em[,sapply(em, as.numeric)]))
