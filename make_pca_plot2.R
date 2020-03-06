#do the pca of ALL
library(data.table)
library(ggplot2)

mesapcs <- fread(file="Z:/data/lauren_mesa/allpops_dosages_joined/mesa_unmerged.eigenvec", header=T)

#pca for the filtered pcs
mesapcs2 <- read.table(file="Z:/data/lauren_mesa/allpops_dosages_joined/ld_pruned_mesa_unmerged.eigenvec", header=T)

ggplot() + geom_point(data=mesapcs,aes(x=PC1,y=PC2)) + theme_bw() + ggtitle("PC1 vs PC2")

ggplot() + geom_point(data=mesapcs2,aes(x=PC1,y=PC2)) + theme_bw() + ggtitle("PC1 vs PC2 geno maf ld filtered")

#check the hapmap pca
hapmap <- fread(file="Z:/data/lauren_mesa/allpops_dosages_joined/hapmap19/hapmap.eigenvec", header=T)

ggplot() + geom_point(data=hapmap,aes(x=PC1,y=PC2)) + theme_bw() + ggtitle("PC1 vs PC2")



