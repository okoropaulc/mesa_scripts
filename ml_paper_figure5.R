#Figure 5

#find the direction of effect of the gene CETP
library(data.table)
library(dplyr)
library(ggplot2)
#read in gencode just to find the gene Id
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)

#CETP gene_id = ENSG00000087237.6

pxcan_adj <- fread(file="Z:/data/twas_mesa/en_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
rf_adj <- fread(file="Z:/data/twas_mesa/rf_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
svr_adj <- fread(file="Z:/data/twas_mesa/svr_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
knn_adj <- fread(file="Z:/data/twas_mesa/knn_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
hdl <- fread(file="Z:/data/twas_mesa/hdl_notrain_rank.txt", header=T)
hdl[,2] <- NULL
names(hdl) <- c("IID", "phenotype")

hdl_en <- inner_join(hdl, pxcan_adj, by = c("IID"="IID"))
en_cetp <- hdl_en[["ENSG00000087237"]]
en_pheno <- hdl_en[["phenotype"]]
en <- data.frame(pheno=en_pheno,cetp=en_cetp)

hdl_rf <- inner_join(hdl, rf_adj, by = c("IID"="IID"))
rf_cetp <- hdl_rf[["ENSG00000087237"]]
rf_pheno <- hdl_rf[["phenotype"]]
rf <- data.frame(pheno=rf_pheno,cetp=rf_cetp)

hdl_svr <- inner_join(hdl, svr_adj, by = c("IID"="IID"))
svr_cetp <- hdl_svr[["ENSG00000087237"]]
svr_pheno <- hdl_svr[["phenotype"]]
svr <- data.frame(pheno=svr_pheno,cetp=svr_cetp)


#plot them
#CETP
library(ggplot2)
ggplot(en, aes(x=cetp, y=pheno)) + geom_point(lwd=3)  +
  xlab("predicted gene expression") + ylab("HDL (rank normalized)") + #xlim(-1.15,0.15) +
  theme_classic(20) + ggtitle("Elastic Net Gene CETP") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

p1 <- ggplot(en, aes(x=cetp, y=pheno)) + geom_point(lwd=3)  +
  xlab("predicted expression") + ylab("HDL (rank normalized)") +
  theme_classic(50) + ggtitle("Elastic Net") + geom_smooth(method="lm", color="red", lwd=3, se=TRUE) +
  geom_density_2d(lwd=2.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")


p2 <- ggplot(rf, aes(x=cetp, y=pheno)) + geom_point(lwd=3)  +
  xlab("predicted expression") + ylab("HDL (rank normalized)") +
  theme_classic(50) + ggtitle("Random Forest") + geom_smooth(method="lm", color="red", lwd=3, se=TRUE) +
  geom_density_2d(lwd=2.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")


p3 <- ggplot(svr, aes(x=cetp, y=pheno)) + geom_point(lwd=3)  +
  xlab("predicted expression") + ylab("HDL (rank normalized)") +
  theme_classic(50) + ggtitle("Support Vector") + geom_smooth(method="lm", color="red", lwd=3, se=TRUE) +
  geom_density_2d(lwd=2.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

library(gridExtra)
grid.arrange(p1,p2,p3,nrow=1)

#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
