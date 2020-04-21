library(data.table)
library(dplyr)
"%&%" <- function(a,b) paste(a,b, sep = "")
library(ggplot2)

gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)

#cetp <- "ENSG00000087237"
#st8sia4 <- "ENSG00000113532"

pop <- "" #this means ALL
trait <- "hdl"

#for (trait in c("hdl", "ldl", "chol", "trig", "plt5")){}

#for (pop in c("AFA", "CAU", "HIS")){}

en <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "en_rank_" %&% trait %&% "_assoc.txt", header=T)
rf <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "rf_rank_" %&% trait %&% "_assoc.txt", header=T)
svr <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "svr_rank_" %&% trait %&% "_assoc.txt", header=T)
knn <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "knn_rank_" %&% trait %&% "_assoc.txt", header=T)


#prepare data for manhattan plot

gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

gencode1 <- gencode[,c(1,2,3,4)]

en_man <- inner_join(gencode1, en, by = c("gene_id" = "gene"))
en_man <- en_man[,c(1,3,4,7)]#drop gene_id and use gene_name to help annotation in the plot

rf_man <- inner_join(gencode1, rf, by = c("gene_id" = "gene"))
rf_man <- rf_man[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svr_man <- inner_join(gencode1, svr, by = c("gene_id" = "gene"))
svr_man <- svr_man[,c(1,3,4,7)] #drop gene_id

knn_man <- inner_join(gencode1, knn, by = c("gene_id" = "gene"))
knn_man <- knn_man[,c(1,3,4,7)] #drop gene_id


#arrange the data by chrom number and start pos

#take each chrom out, sort it by "start", and rbind them
enassocman <- NULL
rfassocman <- NULL
svrassocman <- NULL
knnassocman <- NULL

#en
for (i in 1:22){
  i <- as.character(i)
  a <- subset(en_man, chr==i)
  a <- a[order(a$start),]
  enassocman <- rbind(enassocman,a)
}
names(enassocman) <- c("CHR", "SNP", "BP", "P")
enassocman$CHR <- as.numeric(enassocman$CHR)
# manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net HDL", suggestiveline = F, ylim=c(0,14),
#           col=c("red", "blue"))



#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
# manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest HDL", suggestiveline = F,
#           col=c("red", "blue"), ylim=c(0,14))



#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
# manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector HDL", suggestiveline = F, ylim=c(0,14),
#           col=c("red", "blue"))


#knn
for (i in 1:22){
  i <- as.character(i)
  a <- subset(knn_man, chr==i)
  a <- a[order(a$start),]
  knnassocman <- rbind(knnassocman,a)
}
names(knnassocman) <- c("CHR", "SNP", "BP", "P")
knnassocman$CHR <- as.numeric(knnassocman$CHR)
#knnassocman$P <- as.numeric(knnassocman$P)
# manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN HDL", suggestiveline = F, ylim=c(0,14),
#           col=c("red", "blue"))

fwrite(enassocman, file="Z:/ml_paper_figs/en_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassocman, file="Z:/ml_paper_figs/rf_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassocman, file="Z:/ml_paper_figs/svr_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassocman, file="Z:/ml_paper_figs/knn_assoc_manplot_df.txt", row.names=F, quote=F, sep="\t")


#####Make the manhattan plot
#this script for mamhattan plot below was gotten from
#https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
#I only made a few changes
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

en_df.tmp <- enassocman %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% 
  left_join(enassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 1e-6, "yes", "no"))

en_axisdf <- en_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

rf_df.tmp <- rfassocman %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% 
  left_join(rfassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 1e-6, "yes", "no"))

rf_axisdf <- rf_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

svr_df.tmp <- svrassocman %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% 
  left_join(svrassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 1e-6, "yes", "no"))

svr_axisdf <- svr_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

knn_df.tmp <- knnassocman %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% 
  left_join(knnassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 1e-6, "yes", "no"))

knn_axisdf <- knn_df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

en_gg <- ggplot(en_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + scale_x_continuous( label = en_axisdf$CHR, breaks= en_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("Elastic Net")) + labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-8), colour="red") +
  #geom_hline(yintercept = -log10(1e-6), linetype="dashed") +
  geom_label_repel(data=en_df.tmp[en_df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) + 
  theme_classic(base_size = 22) + theme(legend.position = "none", axis.text.x = element_blank()) + ylab(expression(~~-log[10](italic(p))))


rf_gg <- ggplot(rf_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + scale_x_continuous( label = rf_axisdf$CHR, breaks= rf_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("Random Forest")) + labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-8), colour="red") +
  #geom_hline(yintercept = -log10(1e-6), linetype="dashed") +
  geom_label_repel(data=rf_df.tmp[rf_df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) + 
  theme_classic(base_size = 22) + theme(legend.position = "none", axis.text.x = element_blank()) + ylab(expression(~~-log[10](italic(p))))

svr_gg <- ggplot(svr_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + scale_x_continuous( label = svr_axisdf$CHR, breaks= svr_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("Support Vector")) + labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-8), colour="red") +
  #geom_hline(yintercept = -log10(1e-6), linetype="dashed") +
  geom_label_repel(data=svr_df.tmp[svr_df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) + 
  theme_classic(base_size = 22) + theme(legend.position = "none", axis.text.x = element_blank()) + ylab(expression(~~-log[10](italic(p))))

knn_gg <- ggplot(knn_df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("red","blue"), 22 )) + scale_x_continuous( label = knn_axisdf$CHR, breaks= knn_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14), breaks = seq(0,14,2)) + ggtitle(paste0("K Nearest Neighbor")) + labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(5e-8), colour="red") +
  #geom_hline(yintercept = -log10(1e-6), linetype="dashed") +
  #geom_label_repel(data=knn_df.tmp[knn_df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) + 
  theme_classic(base_size = 22) + theme(legend.position = "none", axis.text.x = element_blank()) + ylab(expression(~~-log[10](italic(p))))


#plot all 4 plots in one page
library(gridExtra)
#fig <- grid.arrange(en_gg,rf_gg,svr_gg,knn_gg,nrow=2)
#print(fig)
tiff("/Users/okoro/OneDrive/Desktop/mesa_scripts/ml_paper_figs/Fig5.tiff", width = 24, height = 16, units = 'cm', res = 300, compression = 'lzw')
grid.arrange(en_gg,rf_gg,svr_gg,knn_gg,nrow=2)
dev.off()
