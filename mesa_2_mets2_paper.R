#ALL
library(data.table)
library(dplyr)
en_all <- read.table(file="Z:/data/mesa_models/split_mesa/results/all_chr_AFA_model_summaries.txt", header=T)#it has no header
en_all$gene_id <- as.character(en_all$gene_id)
en_all <- subset(en_all, en_all$cv_R2_avg >= 0.01)
en_all <- en_all[,c(1,2,10)]
for (i in 1:length(en_all$gene_id)){
  en_all$gene_id[i] <- gsub('\\.[0-9]+','',en_all$gene_id[i])
} #just to remove the decimal places in the gene_id

en_all_2_mets <- read.table(file = "Z:/data/mesa_models/mets_dosages/new_predixcan_AFA2METS_cor", header = T)
en_all_2_mets$gene <- as.character(en_all_2_mets$gene)
filt_en_all_2_mets <- inner_join(en_all, en_all_2_mets, by = c("gene_id" = "gene")) #filter by CV R2 > 0.01
mean(filt_en_all_2_mets$spearman)
#filt_en_all_2_mets <- subset(filt_en_all_2_mets, spearman >= 0.1)

fwrite(filt_en_all_2_mets, file="Z:/data/ml_paper/en_all_2_mets_rho_all.txt", sep="\t", row.names=F, quote=F)


rf_afhi <- fread(file = "Z:/data/mesa_models/python_ml_models/merged_chunk_results/best_grid_rf_all_chrom.txt", header = T)
rf_afhi$Gene_ID <- as.character(rf_afhi$Gene_ID)
rf_afhi <- subset(rf_afhi, rf_afhi$CV_R2 >= 0.01)
for (i in 1:length(rf_afhi$Gene_ID)){
  rf_afhi$Gene_ID[i] <- gsub('\\.[0-9]+','',rf_afhi$Gene_ID[i])
} #just to remove the decimal places in the gene_id

rf_afhi_2_mets <- read.table(file = "Z:/data/mesa_models/python_ml_models/results/grid_optimized_AFA_2_METS_rf_cor_test_full_chr.txt", header = T)
rf_afhi_2_mets$gene_id <- as.character(rf_afhi_2_mets$gene_id)

filt_rf_afhi_2_mets <- inner_join(rf_afhi, rf_afhi_2_mets, by = c("Gene_ID" = "gene_id")) #filter by CV R2 > 0.01
filt_rf_afhi_2_mets <- filt_rf_afhi_2_mets[,c(1,2,3,13)]
names(filt_rf_afhi_2_mets) <- c("gene", "gene_name", "CV_R2", "spearman")
mean(filt_rf_afhi_2_mets$spearman)
#filt_rf_afhi_2_mets <- subset(filt_rf_afhi_2_mets, spearman >= 0.1)
fwrite(filt_rf_afhi_2_mets, file="Z:/data/ml_paper/rf_all_2_mets_rho_all.txt", sep="\t", row.names=F, quote=F)




#########################################################################################################################################
#compare ALL rho of EN against RF, SVR, KNN

#read in en
en_all <- fread(file="Z:/data/ml_paper/en_all_2_mets_rho_0.1.txt", header=T)
en_all$gene_id <- as.character(en_all$gene_id)
rf_all <- fread(file="Z:/data/ml_paper/rf_all_2_mets_rho_0.1.txt", header=T)
rf_all$gene <- as.character(rf_all$gene)
svr_all <- fread(file="Z:/data/ml_paper/svr_all_2_mets_rho_0.1.txt", header=T)
svr_all$gene <- as.character(svr_all$gene)
knn_all <- fread(file="Z:/data/ml_paper/knn_all_2_mets_rho_0.1.txt", header=T)
knn_all$gene <- as.character(knn_all$gene)

en_rf <- inner_join(en_all, knn_all, by = c("gene_id"="gene"))
en_rf <- en_rf[,c(4,7)]
names(en_rf) <- c("en","rf")
cor.test(en_rf$en, en_rf$rf)
library(ggplot2) #i used geom smooth to add the red regression line
ggplot(en_rf, aes(x=en, y=rf)) + geom_point() + geom_abline(intercept=0, slope=1, color="blue") +
  geom_smooth(method="lm", color="red") + xlab("Elastic Net") + ylab("K-Nearest Neighbor") +ylim(0.1,0.85) + xlim(0.1,.85)+ #(-0.1,1) + ylim(-0.1,1) +
  theme_classic(20)# + ggtitle("Cross Validation Perfromance") (save at width=850, height=600) + xlim(-1,1) + ylim(-1,1)



#find the unique genes in each algorithm
enrf_en <- anti_join(en_all, rf_all, by = c("gene_id"="gene")) #1092 genes
enrf_rf <- anti_join(rf_all, en_all, by = c("gene"="gene_id")) #292 genes

ensvr_en <- anti_join(en_all, svr_all, by = c("gene_id"="gene")) #1155 genes
ensvr_svr <- anti_join(svr_all, en_all, by = c("gene"="gene_id")) #296 genes

enknn_en <- anti_join(en_all, knn_all, by = c("gene_id"="gene")) #1564 genes
enknn_knn <- anti_join(knn_all, en_all, by = c("gene"="gene_id")) #211 genes

#make bar plot
comp_bar <- data.frame(comb=c("EN_vs_RF","EN_vs_RF","EN_vs_SVR",
                              "EN_vs_SVR","EN_vs_KNN","EN_vs_KNN"), 
                       alg=c("enrf_en","enrf_rf","ensvr_en","ensvr_svr","enknn_en","enknn_knn"),
                       count=c(1092,292,1155,296,1564,211))

+ geom_text(aes(label=no), vjust=3, color="white", size=10) +
  scale_fill_manual(values = c("blue","violet","red")) + 
  theme_minimal(20) # scale_fill_manual is how to give my own color specification

library(ggplot2)
ggplot(data=comp_bar, aes(x=comb, y=count, fill=alg)) + #scale_fill_manual(values = c("red","violet","red", "green", "red","blue")) + 
  geom_bar(stat="identity", position=position_dodge()) + theme_minimal(20) + 
  geom_text(aes(label=count),vjust=2,color="white", position = position_dodge(0.9), size=6) + xlab("Pairs")

## Make violin plot for all
en_all <- fread(file="Z:/data/ml_paper/en_all_2_mets_rho_all.txt", header=T)
en_all$gene_id <- as.character(en_all$gene_id)
rf_all <- fread(file="Z:/data/ml_paper/rf_all_2_mets_rho_all.txt", header=T)
rf_all$gene <- as.character(rf_all$gene)
svr_all <- fread(file="Z:/data/ml_paper/svr_all_2_mets_rho_all.txt", header=T)
svr_all$gene <- as.character(svr_all$gene)
knn_all <- fread(file="Z:/data/ml_paper/knn_all_2_mets_rho_all.txt", header=T)
knn_all$gene <- as.character(knn_all$gene)


en_all_2mets <- data.frame(spearman=en_all$spearman, 
                           Algorithm=rep("EN", length(en_all$spearman)))

rf_all_2mets <- data.frame(spearman=rf_all$spearman, 
                           Algorithm=rep("RF", length(rf_all$spearman)))

svr_all_2mets <- data.frame(spearman=svr_all$spearman, 
                            Algorithm=rep("SVR", length(svr_all$spearman)))

knn_all_2mets <- data.frame(spearman=knn_all$spearman, 
                            Algorithm=rep("KNN", length(knn_all$spearman)))

alg_mesa_2_mets <- rbind(en_all_2mets, rf_all_2mets, svr_all_2mets, knn_all_2mets)

library(ggplot2)
ggplot(alg_mesa_2_mets, aes(x=Algorithm, y=spearman, color=Algorithm, fill=Algorithm)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.2, color="black") +
  theme_bw(20) #+ ylim(-1,1)


ggplot(alg_mesa_2_mets, aes(x=Algorithm, y=spearman, color=Algorithm, fill=Algorithm)) + 
  geom_violin(trim=F, alpha=0.6) + geom_boxplot(width=0.2, color="black", lwd=1.2) +
  stat_summary(fun.y=mean, geom = "point", size=3, color="white") + theme_linedraw(20) +
  scale_y_continuous(breaks=seq(-1.0, 1.0, 0.2), limits=c(-1.0, 1.0)) +
  theme(panel.grid.minor = element_blank()) + xlab("Elastic Net")
