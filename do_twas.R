#
#rank normalize the lipid trait and also do twas
library(data.table)
library(dplyr)
"%&%" <- function(a,b) paste(a,b, sep = "")

#hdl <- fread(file="Z:/data/twas_mesa/hdl_notrain.txt", header=T)
#ldl <- fread(file="Z:/data/twas_mesa/ldl_notrain.txt", header=T)
#chol <- fread(file="Z:/data/twas_mesa/chol_notrain.txt", header=T)
#trig <- fread(file="Z:/data/twas_mesa/trig_notrain.txt", header=T)

#rankinv <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

#hdl$rank_hdl <- rankinv(hdl$hdl)
#ldl$rank_ldl <- rankinv(ldl$ldl)
#chol$rank_chol <- rankinv(chol$chol)
#trig$rank_trig <- rankinv(trig$trig)

#write.table(trig, file="Z:/data/twas_mesa/trig_notrain_rank.txt",quote=F, row.names=F, sep="\t")
#write.table(ldl, file="Z:/data/twas_mesa/ldl_notrain_rank.txt",quote=F, row.names=F, sep="\t")
#write.table(hdl, file="Z:/data/twas_mesa/hdl_notrain_rank.txt",quote=F, row.names=F, sep="\t")
#write.table(chol, file="Z:/data/twas_mesa/chol_notrain_rank.txt",quote=F, row.names=F, sep="\t")


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


# #Read in the rank normalized
# 
# hdl <- fread(file="Z:/data/twas_mesa/hdl_notrain_rank.txt", header=T) #%>% select(c("sidno","rank_" %&% trait))
# ldl <- fread(file="Z:/data/twas_mesa/ldl_notrain_rank.txt", header=T)
# chol <- fread(file="Z:/data/twas_mesa/chol_notrain_rank.txt", header=T)
# trig <- fread(file="Z:/data/twas_mesa/trig_notrain_rank.txt", header=T)
# 
# pop <- "AFA"
# 
#do the twas here
#but first, read in the adjisted expression

for (trait in c("hdl", "ldl", "chol", "trig", "plt5")){
  
  #read in the rank normalized pheno 
  pheno <- fread(file="/home/pokoro/data/twas_mesa/" %&% trait %&% "_notrain_rank.txt", header=T) %>% select(c("sidno","rank_" %&% trait))
  names(pheno) <- c("IID","phenotype")
  
  for (pop in c("AFA", "CAU", "HIS")){
    
    
    #read in the adjusted expressions
    pxcan_adj <- fread(file="/home/pokoro/data/twas_mesa/en_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
    pxcangenes <- colnames(pxcan_adj)[2:length(pxcan_adj)]
    rf_adj <- fread(file="/home/pokoro/data/twas_mesa/rf_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
    rfgenes <- colnames(rf_adj)[2:length(rf_adj)]
    svr_adj <- fread(file="/home/pokoro/data/twas_mesa/svr_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
    svrgenes <- colnames(svr_adj)[2:length(svr_adj)]
    knn_adj <- fread(file="/home/pokoro/data/twas_mesa/knn_" %&% pop %&% "_allchrom_r2_0.01_adj_pred_exp.txt", header = T)
    knngenes <- colnames(knn_adj)[2:length(knn_adj)]
    
    
    #reading the phenotypes with no training samples
    enpheno <- inner_join(pheno, pxcan_adj, by = c("IID"="IID"))
    enpheno$IID <- NULL
    rfpheno <- inner_join(pheno, rf_adj, by = c("IID"="IID"))
    rfpheno$IID <- NULL
    svrpheno <- inner_join(pheno, svr_adj, by = c("IID"="IID"))
    svrpheno$IID <- NULL
    knnpheno <- inner_join(pheno, knn_adj, by = c("IID"="IID"))
    knnpheno$IID <- NULL
    
    #do twas
    pxcanassoc <- association(enpheno, pxcangenes)
    #pxcanassoc$p <- as.numeric(pxcanassoc$p)
    rfassoc <- association(rfpheno, rfgenes)
    #rfassoc$p <- as.numeric(rfassoc$p)
    svrassoc <- association(svrpheno, svrgenes)
    #svrassoc$p <- as.numeric(svrassoc$p)
    knnassoc <- association(knnpheno, knngenes)
    #knnassoc$p <- as.numeric(knnassoc$p)
    
    #write out the plt assoc
    fwrite(pxcanassoc, file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_en_rank_" %&% trait %&% "_assoc.txt", row.names=F, quote=F, sep="\t")
    fwrite(rfassoc, file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_rf_rank_" %&% trait %&% "_assoc.txt", row.names=F, quote=F, sep="\t")
    fwrite(svrassoc, file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_svr_rank_" %&% trait %&% "_assoc.txt", row.names=F, quote=F, sep="\t")
    fwrite(knnassoc, file="/home/pokoro/data/twas_mesa/" %&% pop %&% "_knn_rank_" %&% trait %&% "_assoc.txt", row.names=F, quote=F, sep="\t")
    
    
    
  }
  
}

