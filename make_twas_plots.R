library(data.table)
library(dplyr)
"%&%" <- function(a,b) paste(a,b, sep = "")
library(ggplot2)

gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)

#cetp <- "ENSG00000087237"
#st8sia4 <- "ENSG00000113532"

pop <- "HIS_"
trait <- "hdl"

#for (trait in c("hdl", "ldl", "chol", "trig", "plt5")){}

#for (pop in c("AFA", "CAU", "HIS")){}

en <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "en_rank_" %&% trait %&% "_assoc.txt", header=T)
rf <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "rf_rank_" %&% trait %&% "_assoc.txt", header=T)
svr <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "svr_rank_" %&% trait %&% "_assoc.txt", header=T)
knn <- fread(file="Z:/data/twas_mesa/" %&% pop %&% "knn_rank_" %&% trait %&% "_assoc.txt", header=T)


#Make manhattan plot for the gene pheno assoc results
#library(qqman)
#see example table for manhattan plot
#eggwas <- gwasResults #this comes with qqman

#gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

#gencode1 <- subset(gencode, gene_type=="protein_coding")
#gencode1$chr <- as.numeric(gencode1$chr)
#chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
gencode1 <- gencode[,c(1,2,3,4)]

en_man <- inner_join(gencode1, en, by = c("gene_id" = "gene"))
en_man <- en_man[,c(1,3,4,7)]#drop gene_id and use gene_name to help annotation in the plot

rf_man <- inner_join(gencode1, rf, by = c("gene_id" = "gene"))
rf_man <- rf_man[,c(1,3,4,7)] #drop gene_id and use gene_name to help annotation in the plot

svr_man <- inner_join(gencode1, svr, by = c("gene_id" = "gene"))
svr_man <- svr_man[,c(1,3,4,7)] #drop gene_id

knn_man <- inner_join(gencode1, knn, by = c("gene_id" = "gene"))
knn_man <- knn_man[,c(1,3,4,7)] #drop gene_id


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




#Make qq plot with the P-values
p <- enassocman[,4]
logp<- -log10(p)
logp <- sort(logp, decreasing = T, na.last = NA) %>% as.data.frame()

create_qq_input<-function(df, limit=nrow(df),range=nrow(df)){ #takes in one column df of pvalues - assumes no NAs, can also limit to the top n hits. Default of this will be just all pvalues
  print("Size of parental set is presumed to be " %&% range)
  print("Provided number of pvalues is " %&% nrow(df))
  print("limiting output to top " %&% limit %&% " rows")
  count <- nrow(df) #count number of pvals that are not na
  ExpP <- -log10((1:count)/(range+1)) %>% as.data.frame() #generate a list of expected pvalues equal in length to the obs pvalues that are not na
  qqvals<-cbind.data.frame(df[1:limit,],ExpP[1:limit,])
  colnames(qqvals)<-c("Observed","Expected")
  return(qqvals)
}

qq_input<-create_qq_input(logp)

qq1 <- ggplot(qq_input,aes(x=Expected,y=Observed)) + 
  geom_point(shape=1) + 
  #coord_cartesian(xlim=c(-0.05,8.05),ylim=c(-0.05,30.05)) + 
  theme_bw(12) + 
  scale_colour_manual(values=c(rgb(163,0,66,maxColorValue = 255),"dark gray")) + 
  theme(legend.position = c(0.01,0.99),legend.justification = c(0,1)) + 
  geom_abline(slope=1,intercept = 0) + 
  xlab(expression(Expected~~-log[10](italic(p)))) + 
  ylab(expression(Observed~~-log[10](italic(p)))) +
  expand_limits(x = 0, y = 0) + ylim(0,14)




#Make manhattan plot with the P-Values
enassocman$CHR <- as.factor(enassocman$CHR)
enassocman$logp <- -log10(enassocman$P)

ggplot(enassocman, aes(x=CHR,y=logp)) + geom_point()



#https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
#Link where this code was copied from
# Variables ====
mypalette <- c("red", "blue")#c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
#mysnps <- c("rs11801961","rs116558464","rs61703161") # snps to highlight
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims, breaks = seq(ylims[1],ylims[2],2)) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig), colour="red") +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_classic(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}


gg.manhattan(enassocman, threshold=1e-6, hlight=NA, col=mypalette, ylims=c(0,14), title="EN")
gg.manhattan(rfassocman, threshold=1e-6, hlight=NA, col=mypalette, ylims=c(0,14), title="RF")
gg.manhattan(svrassocman, threshold=1e-6, hlight=NA, col=mypalette, ylims=c(0,14), title="SVR")
gg.manhattan(knnassocman, threshold=1e-6, hlight=NA, col=mypalette, ylims=c(0,14), title="KNN")




df.tmp <- knnassocman %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% 
  left_join(knnassocman, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot) %>%
  mutate( is_highlight=ifelse(SNP %in% NA, "yes", "no")) %>% mutate( is_annotate=ifelse(P < 1e-6, "yes", "no"))

data=df.tmp[df.tmp$is_annotate=="yes",]

axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g1 <-ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(mypalette, 22 )) + scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,14)) + ggtitle(paste0("EN")) + labs(x = "Chromosome") +
  geom_hline(yintercept = -log10(sig), colour="red") #+
  #geom_hline(yintercept = -log10(sugg), linetype="dashed") #+
  #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3)

g2 <- ifelse(nrow(data) > 0, print(g1+geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], 
                                                       aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3)),
             print(g1)) + theme_classic(base_size = 22) + theme(plot.title = element_text(hjust = 0.5),legend.position="none",
                                                                panel.border = element_blank(),
                                                                panel.grid.major.x = element_blank(),
                                                                panel.grid.minor.x = element_blank())

g1 + theme_classic(base_size = 22) + theme(plot.title = element_text(hjust = 0.5),legend.position="none",
                                           panel.border = element_blank(),
                                           panel.grid.major.x = element_blank(),
                                           panel.grid.minor.x = element_blank())


#plot all 4 plots in one page
library(gridExtra)
grid.arrange(p1,p2,p3,p4,nrow=2)