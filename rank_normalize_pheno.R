#
#rank normalize the lipid trait and also do twas
library(data.table)
library(dplyr)

hdl <- fread(file="Z:/data/twas_mesa/hdl_notrain.txt", header=T)
ldl <- fread(file="Z:/data/twas_mesa/ldl_notrain.txt", header=T)
chol <- fread(file="Z:/data/twas_mesa/chol_notrain.txt", header=T)
trig <- fread(file="Z:/data/twas_mesa/trig_notrain.txt", header=T)

rankinv <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

hdl$rank_hdl <- rankinv(hdl$hdl)
ldl$rank_ldl <- rankinv(ldl$ldl)
chol$rank_chol <- rankinv(chol$chol)
trig$rank_trig <- rankinv(trig$trig)

write.table(trig, file="Z:/data/twas_mesa/trig_notrain_rank.txt",quote=F, row.names=F, sep="\t")
write.table(ldl, file="Z:/data/twas_mesa/ldl_notrain_rank.txt",quote=F, row.names=F, sep="\t")
write.table(hdl, file="Z:/data/twas_mesa/hdl_notrain_rank.txt",quote=F, row.names=F, sep="\t")
write.table(chol, file="Z:/data/twas_mesa/chol_notrain_rank.txt",quote=F, row.names=F, sep="\t")

#do the twas here
#but first, read in the adjisted expression

#writeout the adjusted expressions
pxcan_adj <- fread(file="Z:/data/twas_mesa/en_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
rf_adj <- fread(file="Z:/data/twas_mesa/rf_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
svr_adj <- fread(file="Z:/data/twas_mesa/svr_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
knn_adj <- fread(file="Z:/data/twas_mesa/knn_ALL_allchrom_r2_0.01_adj_pred_exp.txt", header = T)



#ascociation function

association <- function(merged, genes, test_type = "linear") {
  assoc_df <- NULL # Init association dataframe
  
  # Perform test between each pred_gene_exp column and phenotype
  for (gene in genes) {
    pred_gene_exp <- merged[[gene]]
    if (test_type == "logistic") { 
      model <- glm(phenotype ~ pred_gene_exp, data = merged, family = binomial)
    } else if (test_type == "linear") {
      model <- lm(phenotype ~ pred_gene_exp, data = merged)
    } else if (test_type == "survival") {
      # TODO: survival analysis
      model <- NULL
    }
    results <- coef(summary(model))[c(2,6,8,4)]
    line <- c(gene,results)
    assoc_df <- rbind(assoc_df,line)
  }
  
  # Specify column names of assoc_df
  if (test_type == "logistic") {
    colnames(assoc_df) <- c("gene", "beta", "z", "p", "se(beta)")
  } else if (test_type == "linear") {
    colnames(assoc_df) <- c("gene", "beta", "t", "p", "se(beta)")
  } else if (test_type == "survival") {
    # TODO
  }
  return(as.data.frame(assoc_df))
}

#reading the phenotypes with no training samples
#platelet count
plt <- fread(file="Z:/data/twas_mesa/rankplt5_notrain.txt", header=T)
names(plt) <- c("IID","phenotype")

enplt <- inner_join(plt, pxcan_adj, by = c("IID"="IID"))
enplt$IID <- NULL
rfplt <- inner_join(plt, rf_adj, by = c("IID"="IID"))
rfplt$IID <- NULL
svrplt <- inner_join(plt, svr_adj, by = c("IID"="IID"))
svrplt$IID <- NULL
knnplt <- inner_join(plt, knn_adj, by = c("IID"="IID"))
knnplt$IID <- NULL

#do twas
pxcanassoc <- association(enplt, pxcangenes)
#pxcanassoc$p <- as.numeric(pxcanassoc$p)
rfassoc <- association(rfplt, rfgenes)
#rfassoc$p <- as.numeric(rfassoc$p)
svrassoc <- association(svrplt, svrgenes)
#svrassoc$p <- as.numeric(svrassoc$p)
knnassoc <- association(knnplt, knngenes)
#knnassoc$p <- as.numeric(knnassoc$p)

#write out the plt assoc
fwrite(pxcanassoc, file="Z:/data/twas_mesa/en_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(rfassoc, file="Z:/data/twas_mesa/rf_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(svrassoc, file="Z:/data/twas_mesa/svr_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")
fwrite(knnassoc, file="Z:/data/twas_mesa/knn_rankplt5_assoc.txt", row.names=F, quote=F, sep="\t")


#don't mind the variable here. I just used same variable for different files.
#read in the assoc again so R can treat the pvalue as numeric
assoc_en <- fread(file="Z:/data/twas_mesa/en_rank_trig_assoc.txt",header=T)
assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rank_trig_assoc.txt",header=T)
assoc_svr <- fread(file="Z:/data/twas_mesa/svr_rank_trig_assoc.txt",header=T)
assoc_knn <- fread(file="Z:/data/twas_mesa/knn_rank_trig_assoc.txt",header=T)

#ggplot after filtering
library(qqman)
qq(assoc_en$p, main="Elastic Net Rank_TRIG")
qq(assoc_rf$p, main="Random Forest Rank_TRIG")
qq(assoc_svr$p, main="Support Vector Rank_TRIG")
qq(assoc_knn$p, main="KNN Rank_TRIG")
#rf_pheno_qq_filtered


#Make manhattan plot for the gene pheno assoc results
library(qqman)
#see example table for manhattan plot
#eggwas <- gwasResults #this comes with qqman

#manhattan for assoc filtered by cv r2
#merge the gene pheno assoc results with gencode annotation
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
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

#en
for (i in 1:22){
  i <- as.character(i)
  a <- subset(en_man, chr==i)
  a <- a[order(a$start),]
  enassocman <- rbind(enassocman,a)
}
names(enassocman) <- c("CHR", "SNP", "BP", "P")
enassocman$CHR <- as.numeric(enassocman$CHR)
manhattan(enassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Elastic Net TRIG", annotatePval = -log10(5e-08),
          col=c("red", "blue"))

#rf
for (i in 1:22){
  i <- as.character(i)
  a <- subset(rf_man, chr==i)
  a <- a[order(a$start),]
  rfassocman <- rbind(rfassocman,a)
}
names(rfassocman) <- c("CHR", "SNP", "BP", "P")
rfassocman$CHR <- as.numeric(rfassocman$CHR)
manhattan(rfassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Random Forest TRIG", annotatePval = -log10(5e-08),
          col=c("red", "blue"))

#svr
for (i in 1:22){
  i <- as.character(i)
  a <- subset(svr_man, chr==i)
  a <- a[order(a$start),]
  svrassocman <- rbind(svrassocman,a)
}
names(svrassocman) <- c("CHR", "SNP", "BP", "P")
svrassocman$CHR <- as.numeric(svrassocman$CHR)
manhattan(svrassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="Support Vector TRIG",annotatePval = -log10(5e-08),
          col=c("red", "blue"))

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
manhattan(knnassocman, chr="CHR", bp="BP", snp="SNP", p="P", main="KNN TRIG", annotatePval=-log10(5e-08),
          col=c("red", "blue"))



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
ggplot(en, aes(x=cetp, y=pheno)) + geom_point()  +
  xlab("predicted gene expression") + ylab("HDL (rank normalized)") + #xlim(-1.15,0.15) +
  theme_classic(20) + ggtitle("Elastic Net Gene CETP") #(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")

ggplot(en, aes(x=cetp, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("HDL (rank normalized)") +
  theme_classic(20) + ggtitle("Elastic Net Gene CETP") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) +
  geom_density_2d(lwd=1.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")


ggplot(rf, aes(x=cetp, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("HDL (rank normalized)") +
  theme_classic(20) + ggtitle("Random Forest Gene CETP") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) +
  geom_density_2d(lwd=1.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")


ggplot(svr, aes(x=cetp, y=pheno)) + geom_point()  +
  xlab("predicted expression") + ylab("HDL (rank normalized)") +
  theme_classic(20) + ggtitle("Support Vector Gene CETP") + geom_smooth(method="lm", color="red", lwd=2, se=TRUE) +
  geom_density_2d(lwd=1.5)#(save at width=850, height=600) + xlim(0,1) + ylim(0,1)
#+ geom_abline(intercept=0, slope=1, color="blue", lwd=2) + geom_smooth(method="lm", color="red")




#compare the t statistic performance of EN and RF
#read in the assoc
assoc_en <- fread(file="Z:/data/twas_mesa/en_rank_hdl_assoc.txt",header=T)
assoc_en <- assoc_en[,c(1,3)] #cols are gene and t
names(assoc_en) <- c("gene","en")

assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rank_hdl_assoc.txt",header=T)
assoc_rf <- assoc_rf[,c(1,3)] #cols are gene and t
names(assoc_rf) <- c("gene","rf")

#join both and keep all rows
en_rf <- full_join(assoc_en, assoc_rf, by = c("gene"="gene"))
#replace all NAs with 0
en_rf[is.na(en_rf)] <- 0

#plot the two t statistic en vs rf
ggplot(en_rf, aes(x=en, y=rf)) + geom_point()  +
  xlab("Elastic Net") + ylab("Random Forest") + #xlim(-1.15,0.15) +
  theme_classic(20) + ggtitle("EN vs RF on HDL") + geom_abline(intercept=0, slope=1, color="blue") #+ 
  #geom_smooth(method="lm", color="red")
#(save at width=900, height=650) + xlim(0,1) + ylim(0,1)

x <- en_rf$en
y <- en_rf$rf

plot(x, y, type="p", xaxt='n', yaxt='n')
arrows(min(x), 0, max(x), 0, lwd=1, length=0.15)
arrows(0, min(y), 0, max(y), lwd=1, length=0.15)

text(0, min(y), "y", pos=2)
text(min(x), 0, "x", pos=3)


#use function from stackexchange
#https://stackoverflow.com/questions/17753101/center-x-and-y-axis-with-ggplot2

theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
                           color = "black", size = 1, 
                           xlab = "x", ylab = "y",
                           ticks = 10,
                           textsize = 3,
                           xlimit = max(abs(xvals),abs(yvals)),
                           ylimit = max(abs(yvals),abs(xvals)),
                           epsilon = max(xlimit,ylimit)/50){
  
  #INPUT:
  #xvals .- Values of x that will be plotted
  #yvals .- Values of y that will be plotted
  #xgeo  .- x intercept value for y axis
  #ygeo  .- y intercept value for x axis
  #color .- Default color for axis
  #size  .- Line size for axis
  #xlab  .- Label for x axis
  #ylab  .- Label for y axis
  #ticks .- Number of ticks to add to plot in each axis
  #textsize .- Size of text for ticks
  #xlimit .- Limit value for x axis 
  #ylimit .- Limit value for y axis
  #epsilon .- Parameter for small space
  
  
  #Create axis 
  xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))
  
  #Add axis
  theme.list <- 
    list(
      theme_void(), #Empty the current theme
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
      annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
      annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
      xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
      ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
    )
  
  #Add ticks programatically
  ticks_x <- round(seq(-xlimit, xlimit, length.out = ticks),2)
  ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)
  
  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:ticks){
    
    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))
    
    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))
    
    #Add ticks to geom line for x axis
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
                                             data = xtick, size = size, 
                                             color = color)
    
    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
                                            x = ticks_x[k], 
                                            y = ygeo - 2.5*epsilon,
                                            size = textsize,
                                            label = paste(ticks_x[k]))
    
    
    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = size, 
                                             color = color)
    
    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                          x = xgeo - 2.5*epsilon, 
                                          y = ticks_y[k],
                                          size = textsize,
                                          label = paste(ticks_y[k]))
  }
  
  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}


ggplot(en_rf) + theme_geometry(en_rf$en, en_rf$rf) +
  geom_point(aes(x = x, y = y), size = 1, color="red") + 
  ggtitle("EN vs RF on HDL") + geom_abline(intercept=0, slope=1, color="blue") + theme_classic(20) + xlab("Elastic Net") + 
  ylab("Random Forest") #save 900 by 650


#color the intersection genes black, and others red
#read in the assoc
assoc_en <- fread(file="Z:/data/twas_mesa/en_rank_hdl_assoc.txt",header=T)
assoc_en <- assoc_en[,c(1,3)] #cols are gene and t
names(assoc_en) <- c("gene","en")

assoc_rf <- fread(file="Z:/data/twas_mesa/rf_rank_hdl_assoc.txt",header=T)
assoc_rf <- assoc_rf[,c(1,3)] #cols are gene and t
names(assoc_rf) <- c("gene","rf")

#join both and keep all rows
en_rf <- full_join(assoc_en, assoc_rf, by = c("gene"="gene"))
en_rf$type <- "" #this type column is so i can keep track of the intercept and unique genes and color it so

#name gene intercept if there is no NA vice versa
for (i in 1:nrow(en_rf)){
  ifelse(is.na(en_rf[i,2]) | is.na(en_rf[i,3]), en_rf$type[i]<-"unique", en_rf$type[i]<-"intersection")
}

#now replace all NAs with 0
en_rf[is.na(en_rf)] <- 0

#plot them
ggplot(en_rf, aes(x=en, y=rf, col=type)) + geom_point()  +
  xlab("Elastic Net") + ylab("Random Forest") + #xlim(-1.15,0.15) +
  theme_classic(20) + ggtitle("EN vs RF on HDL") + geom_abline(intercept=0, slope=1, color="blue") +
  scale_color_manual(values = c("black","red"))#+ 
#geom_smooth(method="lm", color="red")

ggplot(en_rf) + theme_geometry(en_rf$en, en_rf$rf) +
  geom_point(aes(x = x, y = y, col=type), size = 1) + 
  ggtitle("EN vs RF on HDL") + geom_abline(intercept=0, slope=1, color="blue") + theme_classic(20) + xlab("Elastic Net") + 
  ylab("Random Forest") + scale_color_manual(values = c("black","red"))#save 900 by 650




###put gene name

#join both and keep all rows
en_rf <- full_join(assoc_en, assoc_rf, by = c("gene"="gene"))

#put gene name

#merge the gene pheno assoc results with gencode annotation
gencode <- fread(file="Z:/data/mesa_models/gencode.v18.annotation.parsed.txt", header=T)
gencode$gene_id <- as.character(gencode$gene_id)
gencode$gene_type <- as.character(gencode$gene_type)
gencode <- subset(gencode, gene_type=="protein_coding") #this to reduce the df size to only protein coding genes
#remove the decimal in the gene id
for (i in 1:nrow(gencode)){
  gencode$gene_id[i] <- gsub('\\.[0-9]+','',gencode$gene_id[i])
} #just to remove the decimal places in the gene_id

gencode <- gencode[,c(1,2,3)]

en_rf_gene <- inner_join(gencode, en_rf, by = c("gene_id"="gene"))

#replace all NAs with 0
en_rf[is.na(en_rf)] <- 0