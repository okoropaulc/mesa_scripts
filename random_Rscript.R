library(stringr)
hisgex <- read.table(file = "/home/paul/mesa_models/meqtl_sorted_HIS_MESA_Epi_GEX_data_sidno_Nk-10.txt", header = T)
hisgex$PROBE_ID <- as.character(hisgex$PROBE_ID)
for (i in 2:length(hisgex)){
  colnames(hisgex)[i] <- str_sub(names(hisgex[i]),2,6)
}
rfgex <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_HIS_rf_predicted_gene_expr_chr1.txt", header = T, sep = "\t")
rfcol <- colnames(rfgex)[2:352]
rfgex <- rfgex[order(rfgex$X),]
#rfgex <- as.data.frame(t(rfgex))
hisgex <- hisgex[hisgex$PROBE_ID %in% rfcol,]
hisgex <- hisgex[order(hisgex$PROBE_ID),]
#hisgex <- hisgex[,order(names(hisgex)[2:353])]
hisg1 <- subset(hisgex, hisgex$PROBE_ID == "ENSG00000117411.12")
hisg1$PROBE_ID <- NULL
hisg1 <- hisg1[,order(names(hisg1))]

thisg1 <- t(hisg1)
rfg1 <- as.data.frame(rfgex$ENSG00000117411.12)
names(rfg1) <- "rf"

thisg1$rf <- rfgex$ENSG00000117411.12

plot(thisg1, rfg1$rf)

rfgex2 <- read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_HIS_rf_predicted_gene_expr_chr2.txt", header = T, sep = "\t")


"%&%" <- function(a,b) paste(a,b, sep = "")

knn <- NULL
rf <- NULL
svr <- NULL
svrl <- NULL
#el <- NULL

i <- "HIS"

for (chrom in 1:22) {
  no <- as.character(chrom)
  knn <- cbind(knn, read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_knn_predicted_gene_expr_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  rf <- cbind(rf, read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_rf_predicted_gene_expr_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  svr <- cbind(svr, read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_svr_predicted_gene_expr_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  svrl <- cbind(svrl, read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_svrl_predicted_gene_expr_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
  
  #el <- rbind(el, read.table(file = "/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_elnet_cor_test_chr" %&% no %&% ".txt", header = T, stringsAsFactors = F, sep = "\t"))
}

#write out the merged full chromosomes
write.table(knn, file="/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_knn_predicted_gene_expr_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rf, file="/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_rf_predicted_gene_expr_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svr, file="/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_svr_rbf_predicted_gene_expr_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(svrl, file="/home/paul/mesa_models/python_ml_models/results/AFA_2_" %&% i %&% "_svr_linear_predicted_gene_expr_all_chr.txt", quote = FALSE, sep = "\t", row.names = FALSE)


a <- matrix(c(2,1,2,1,2,2,2,2,1,2,1,2), ncol = 3)
a <- as.data.frame(a)
c1 <- read.csv(file = "/home/paul/mesa_models/svr/svr_cis_gt_chr14_CHURC1.csv")
varid <- names(c1)
varid <- varid[1:15]
c2 = c1[, varid, drop=FALSE]
rsid <- c("rs01", "rs02", "rs03", "rs04", "rs05", "rs06", "rs07", "rs08", "rs09", "rs10", "rs11", "rs12", "rs13", "rs14", "rs15")
colnames(c2) <- rsid
gcov <- cov(c2)
