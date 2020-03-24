#plot figure 4 in ML Paper
#manhatan and qqplot
library(data.table)
library(ggplot2)
library(dplyr)
#library(viridis)

#read in the assoc again so R can treat the pvalue as numeric
assoc_en <- fread(file="Z:/data/twas_mesa/en_rank_hdl_assoc.txt",header=T)
assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rank_hdl_assoc.txt",header=T)
assoc_svr <- fread(file="Z:/data/twas_mesa/svr_rank_hdl_assoc.txt",header=T)
assoc_knn <- fread(file="Z:/data/twas_mesa/knn_rank_hdl_assoc.txt",header=T)


library(qqman)
# par("mar") default 5.1,4.1,4.1,2.1
#par("mai") default 1.02,0.82,0.82,0.42
#par("oma") default 0,0,0,0
par(mfrow=c(2,2), mar=c(8.1,7.1,7.1,4.1))#, oma=c(0,0,0,0))#, mai=c(2,2,2,2)) #use this to make all 4 plots into one panel
#par(mar=c(3.1,2.1,2.1,1.1))
qq(assoc_en$p, ylim=c(0,14), cex=3, las=1, cex.lab=3, 
   cex.axis=2, main="Elastic Net", cex.main=4, xlim=c(0,4))

qq(assoc_rf$p, ylim=c(0,14), cex=3, las=1, cex.lab=3, 
   cex.axis=2, main="Random Forest", cex.main=4, xlim=c(0,4))

qq(assoc_svr$p, ylim=c(0,14), cex=3, las=1, cex.lab=3, 
   cex.axis=2, main="Support Vector",cex.main=4, xlim=c(0,4))

qq(assoc_knn$p, ylim=c(0,14), cex=3, las=1, 
   cex.lab=3, cex.axis=2, main="KNN", cex.main=4, xlim=c(0,4))
#main="KNN Rank_HDL", save width=4000, height=3400

#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2)

#merge the gene pheno assoc results with gencode annotation
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

gencode1 <- gencode[,c(1,2,3,4)]

en_man <- inner_join(gencode1, assoc_en, by = c("gene_id" = "gene"))
en_man <- en_man[,c(1,3,4,7)]#drop gene_id and use gene_name to help annotation in the plot

rf_man <- inner_join(gencode1, assoc_rf, by = c("gene_id" = "gene"))
rf_man <- rf_man[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svr_man <- inner_join(gencode1, assoc_svr, by = c("gene_id" = "gene"))
svr_man <- svr_man[,c(1,3,4,7)] #drop gene_id

knn_man <- inner_join(gencode1, assoc_knn, by = c("gene_id" = "gene"))
knn_man <- knn_man[,c(1,3,4,7)] #drop gene_id


#take each chrom out, sort it by "start", and rbind them
enassocman <- NULL
rfassocman <- NULL
svrassocman <- NULL
knnassocman <- NULL

par(mfrow=c(2,2), mar=c(8.1,7.1,7.1,4.1))

#en
for (i in 1:22){
  i <- as.character(i)
  a <- subset(en_man, chr==i)
  a <- a[order(a$start),]
  enassocman <- rbind(enassocman,a)
}
names(enassocman) <- c("CHR", "SNP", "BP", "P")
enassocman$CHR <- as.numeric(enassocman$CHR)
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net", 
          suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=3, cex.axis=1.5, cex.lab=3, cex.main=3)

# cex=2, cex.axis=1.5, cex.lab=4)
#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest", 
          suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=3, cex.axis=1.5, cex.lab=3, cex.main=3)

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector", 
          suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=3, cex.axis=1.5, cex.lab=3, cex.main=3)

#knn
for (i in 1:22){
  i <- as.character(i)
  a <- subset(knn_man, chr==i)
  a <- a[order(a$start),]
  knnassocman <- rbind(knnassocman,a)
}
names(knnassocman) <- c("CHR", "SNP", "BP", "P")
knnassocman$CHR <- as.numeric(knnassocman$CHR)
knnassocman$P <- as.numeric(knnassocman$P)
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN", 
          suggestiveline = F,
          col=c("red", "blue"), ylim=c(0,14), cex=3, cex.axis=1.5, cex.lab=3, cex.main=3)



#Try manhattan plot with GGPLoT
#https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

gwas.dat <- rfassocman
nCHR <- length(unique(gwas.dat$CHR))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(gwas.dat$CHR)){
  nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$BP)
  gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"BP"] + s
  s <- s + nbp[i]
}


axis.set <- gwas.dat %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(gwas.dat$P)))) + 2 
sig <- 5e-8


manhplot <- ggplot(gwas.dat, aes(x = BPcum, y = -log10(P), 
                                 color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

print(manhplot)