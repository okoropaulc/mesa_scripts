#Note:
#That the grid optimized svr was optimized for all possible values of hyperparameter C
#However, the output was named **svr_rbf** especially for the ALL_2_METS
#So ignore the fact that the name contain svr_rbf

library(data.table)
"%&%" <- function(a,b) paste(a,b, sep="")
library(dplyr)


#en_all_2_mets <- fread(file = "Z:/data/mesa_models/mets_dosages/ALL_2_METS_predicted_expression.txt", header = T)

#all2me_rf <- fread(file="Z:/data/mesa_models/python_ml_models/results/grid_optimized_ALL_2_METS_rf_predicted_gene_expr_chr8.txt",header=T)

algs <- c("rf", "svr_rbf", "knn")

for (alg in algs){
  all2mets <- NULL
  for (i in 1:22){
    all2mets <- cbind(all2mets, 
                      fread(file="/home/pokoro/data/mesa_models/python_ml_models/results/grid_optimized_ALL_2_METS_" %&% alg %&% "_predicted_gene_expr_chr" %&% i %&% ".txt",header=T))
    
  }
  fwrite(all2mets, file="/home/pokoro/data/mesa_models/python_ml_models/results/ALL_2_METS_" %&% alg %&% "_predicted_expr_allchrom.txt",
         row.names=F, quote=F, sep="\t")
}


#read in the mbiome diversity index file
diversity <- read.table(file="Z:/paul/mets61/mets61_pheno_alpha_D_index_sorted.txt", header=T)
diversity$Obesity_status <- as.character(diversity$Obesity_status)
#diversity$obese <- ""

for (i in 1:nrow(diversity)){
#  if(diversity$Obesity_status[i]=="Lean"){
#    diversity$obese[i] <- 0
#  }
  ifelse(diversity$Obesity_status[i]=="Lean", diversity$obese[i]<-0, diversity$obese[i]<-1)
}


library(ggplot2)

ggplot(diversity, aes(x=Site, y=Shannon_index, colour=Site)) + geom_violin()

