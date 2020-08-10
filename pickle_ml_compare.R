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



#Check how many genes should be in a model
afachr <- read.table(file="Z:/data/mesa_models/python_ml_models/merged_chunk_results/AFA_best_grid_rf_chr22_full.txt", header=1)
rf <- subset(afachr, CV_R2 > 0.01)











#### do the correlation between pickle predicted and measured expression
#AFA 2 METS with RF
library(dplyr)

chr22 <- read.table("Z:/data/ml_software/METS_RF_predicted_expression_chrom_22_transposed.txt", header=T, stringsAsFactors=F)

for (i in 1:length(chr22)){
  colnames(chr22)[i] <- gsub('\\.[0-9]+','',colnames(chr22)[i])
} #just to remove the decimal places in the gene_id ie columns

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding')){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}
gene_annotation_path = "Z:/data/METS_model/hg19/gencode.v28_annotation.parsed.txt"
gene_anot <- get_gene_annotation(gene_annotation_path, chrom=22)

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}

mets_gex_chr22 <- get_gene_expression("Z:/data/METS_model/hg19/METS_peer10_all_chr_prot_coding_gex.txt", gene_anot)

for (i in 1:length(mets_gex_chr22)){
  colnames(mets_gex_chr22)[i] <- gsub('\\.[0-9]+','',colnames(mets_gex_chr22)[i])
} #just to remove the decimal places in the gene_id ie columns

#select only genes in the predicted expression file
t_chr22 <- as.data.frame(t(chr22))
probe_id <- as.character(rownames(t_chr22))
t_chr22 <- cbind(probe_id, t_chr22)
t_chr22$probe_id <- as.character(t_chr22$probe_id)

t_mets_gex_chr22 <- as.data.frame(t(mets_gex_chr22))
probe_id <- rownames(t_mets_gex_chr22)
t_mets_gex_chr22 <- cbind(probe_id, t_mets_gex_chr22)
t_mets_gex_chr22$probe_id <- as.character(t_mets_gex_chr22$probe_id)

#arrange the columns to be the same
t_chr22 <- t_chr22 %>% select(colnames(t_mets_gex_chr22))

#merge the probe_id
probe_id <- data.frame(probe_id=rownames(t_chr22))
probe_id$probe_id <- as.character(probe_id$probe_id)

t_mets_gex_chr22 <- inner_join(probe_id, t_mets_gex_chr22, by =c("probe_id"="probe_id"))
probe_id <- data.frame(probe_id=t_mets_gex_chr22$probe_id) #just to make sure that t_chr22 has same number of rows in same order
probe_id$probe_id <- as.character(probe_id$probe_id)
t_chr22 <- inner_join(probe_id, t_chr22, by =c("probe_id"="probe_id"))

pred <- as.data.frame(t(t_chr22), stringsAsFactors = F)
colnames(pred) <- t_chr22$probe_id
pred <- pred[-1,] #remove the first row. it is probe_id
obs <- as.data.frame(t(t_mets_gex_chr22), stringsAsFactors = F)
colnames(obs) <- t_mets_gex_chr22$probe_id
obs <- obs[-1,]

#do the spearman correlation between predicted and measured
cor_df <- data.frame(rho=numeric())

for (i in 1:ncol(pred)){
  #since the nrow and ncol are equal and in same order. genes are colums, samples are rows
  test <- cor.test(as.numeric(pred[,i]), as.numeric(obs[,i]), method = "spearman")
  cor_df[i,1] <- as.numeric(test$estimate)
}
cor_df <- cbind(probe_id,cor_df)

#compare with non_pickle cor of chr22
library(data.table)
rf_22 <- fread(file="Z:/data/mesa_models/python_ml_models/results/AFA_2_METS_rf_cor_test_chr22.txt", header=T)
rf_22 <- rf_22[,c(1,9)]
names(rf_22) <- c("probe_id", "non_pickle")

df <- inner_join(cor_df, rf_22, by = c("probe_id"="probe_id"))

#plot it
library(ggplot2)
library(ggpubr)

tiff("/Users/okoro/OneDrive/Desktop/afa_2_mets_rf_chr22_pickle_compare.tiff", width = 16, height = 14, units = "cm",
     res = 300, compression = "lzw")
ggplot(df, aes(x=rho,y=non_pickle)) + geom_point() + xlab("pickled") + ylab("non_pickled") + theme_classic(20) +
  ggtitle("Pickle vs Non_picke in METS") + geom_smooth(method = "lm")
dev.off()

#retain probe_id that's in t_chr22
#t_mets_gex_chr22 <- subset(t_mets_gex_chr22, probe_id %in% t_chr22$probe_id)

# #gencode <- read.table(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
# mets_gex <- read.table(file="Z:/data/METS_model/hg19/METS_peer10_all_chr_prot_coding_gex.txt", header=T)
# mets_gex2 <- mets_gex
# 
# for (i in 1:nrow(mets_gex)){
#   mets_gex$PROBE_ID[i] <- gsub('\\.[0-9]+','',mets_gex$PROBE_ID[i])
# } #just to remove the decimal places in the gene_id ie columns
