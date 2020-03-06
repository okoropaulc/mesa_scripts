#join the populations (AFA, CAU, HIS) genotype into one on their common rsid
#library(dplyr)

"%&%" <- function(a,b) paste(a,b, sep = "")

for (i in 1:22){
  his <- read.table(file="/home/pokoro/data/lauren_mesa/his_dosages/chr" %&% i %&% ".maf0.01.r20.8noambig.dosage.txt.gz", header=F)
  his$V6 <- NULL #aa_freq
  his$V2 <- as.character(his$V2) #rsid
  afa <- read.table(file="/home/pokoro/data/lauren_mesa/afa_dosages/chr" %&% i %&% ".maf0.01.r20.8noambig.dosage.txt.gz", header=F)
  afa$V6 <- NULL
  afa$V2 <- as.character(afa$V2)
  cau <- read.table(file="/home/pokoro/data/lauren_mesa/cau_dosages/chr" %&% i %&% ".maf0.01.r20.8noambig.dosage.txt.gz", header=F)
  cau$V6 <- NULL
  cau$V2 <- as.character(cau$V2)
  
  #column V6 is the aa_freq
  #cols are V1=chr, V2=rsid, V3=pos, V4=ref, V5=alt
  
  
  #find the overlapping snps
  
  his_snps <- his$V2
  afa_snps <- afa$V2
  cau_snps <- cau$V2
  
  snp_join <- c(his_snps, afa_snps, cau_snps)
  snp_join <- unique(snp_join)
  
  snp_join <- intersect(his_snps, afa_snps)
  snp_join <- intersect(snp_join, cau_snps)
  
  his <- subset(his, V2 %in% snp_join)
  afa <- subset(afa, V2 %in% snp_join)
  cau <- subset(cau, V2 %in% snp_join)
  #merge them in this order afa, cau, his with cbind
  #first remove columns V1,V2,V3,V4,V5 from cau and his, since they are now same with the one in afa after taking their snp overlaps
  
  cau[,c(1:5)] <- NULL
  his[,c(1:5)] <- NULL
  
  #join them with cbind afa cau his. Therefore rbind there individual sample.txt files in same order to maintain sample names
  allpops <- cbind(afa,cau,his)
  
  #write out the allpops
  #first on per chromosome
  write.table(allpops, file="/home/pokoro/data/lauren_mesa/allpops_dosages_joined/chr" %&% i %&% ".txt",row.names=F,col.names=F,quote=F,sep="\t")
  #also into whole chrom by using the append=T
  write.table(allpops, file="/home/pokoro/data/lauren_mesa/allpops_dosages_joined/whole_chrom.txt",row.names=F,col.names=F,quote=F,sep="\t",append=T)
  
}

###########create sample.txt file for the new joined pops

#merge the samples.txt in the same order
# #Read in AFA
afa_sam <- read.table(file="/home/pokoro/data/lauren_mesa/afa_dosages/samples.txt", header=F)
# 
# #Read in CAU
cau_sam <- read.table(file="/home/pokoro/data/lauren_mesa/cau_dosages/samples.txt", header=F)
# 
# #Read in CAU
his_sam <- read.table(file="/home/pokoro/data/lauren_mesa/his_dosages/samples.txt", header=F)
# 
allpop_sam <- rbind(afa_sam,cau_sam,his_sam)
# #write out the allpop_sample
write.table(allpop_sam, file="/home/pokoro/data/lauren_mesa/allpops_dosages_joined/samples.txt",row.names=F,col.names=F,quote=F,sep="\t")
# 
# 
# #make spoofed.fam for plink2
# #requires 0 sampleid 0 0 0 0
# 
spoofed <- data.frame(id1=rep(0,nrow(allpop_sam)), id2=rep(0,nrow(allpop_sam)), id3=rep(0,nrow(allpop_sam)),
                       id2=rep(0,nrow(allpop_sam)))
 spoofed <- cbind(allpop_sam, spoofed)
write.table(spoofed, file="/home/pokoro/data/lauren_mesa/allpops_dosages_joined/spoofed.fam", row.names=F,col.names=F,quote=F,sep="\t")















# his <- read.table(file="Z:/data/lauren_mesa/his_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=50)
# his$V6 <- NULL #aa_freq
# his$V2 <- as.character(his$V2) #rsid
# afa <- read.table(file="Z:/data/lauren_mesa/afa_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=50)
# afa$V6 <- NULL
# afa$V2 <- as.character(afa$V2)
# cau <- read.table(file="Z:/data/lauren_mesa/cau_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=50)
# cau$V6 <- NULL
# cau$V2 <- as.character(cau$V2)
# #cau2 <- read.table(file="Z:/data/lauren_mesa/cau_dosages/chr21.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=5)
# 
# #column V6 is the aa_freq
# #cols are V1=chr, V2=rsid, V3=pos, V4=ref, V5=alt
# 
# 
# #find the overlapping snps
# library(dplyr)
# 
# his_snps <- his$V2
# afa_snps <- afa$V2
# cau_snps <- cau$V2
# 
# snp_join <- c(his_snps, afa_snps, cau_snps)
# snp_join <- unique(snp_join)
# 
# snp_join <- intersect(his_snps, afa_snps)
# snp_join <- intersect(snp_join, cau_snps)
# 
# his <- subset(his, V2 %in% snp_join)
# afa <- subset(afa, V2 %in% snp_join)
# cau <- subset(cau, V2 %in% snp_join)
# #merge them in this order afa, cau, his with cbind
# #first remove columns V1,V2,V3,V4,V5 from cau and his, since they are now same with the one in afa after taking their snp overlaps
# 
# cau[,c(1:5)] <- NULL
# his[,c(1:5)] <- NULL
# 
# #join them with cbind afa cau his. Therefore rbind there individual sample.txt files in same order to maintain sample names
# allpops <- cbind(afa,cau,his)
# 
# #write out the allpops
# write.table(allpops, file="Z:/allpops.txt",row.names=F,col.names=F,quote=F,sep="\t")
# 
# 
# #merge the samples.txt in the same order
# #Read in AFA
# #afa_dos <- read.table(file="Z:/data/lauren_mesa/afa_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=3)
# afa_sam <- read.table(file="Z:/data/lauren_mesa/afa_dosages/samples.txt", header=F)
# 
# #Read in CAU
# #cau_dos <- read.table(file="Z:/data/lauren_mesa/cau_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=3)
# cau_sam <- read.table(file="Z:/data/lauren_mesa/cau_dosages/samples.txt", header=F)
# 
# #Read in CAU
# #his_dos <- read.table(file="Z:/data/lauren_mesa/his_dosages/chr22.maf0.01.r20.8noambig.dosage.txt.gz", header=F, nrows=3)
# his_sam <- read.table(file="Z:/data/lauren_mesa/his_dosages/samples.txt", header=F)
# 
# allpop_sam <- rbind(afa_sam,cau_sam,his_sam)
# #write out the allpop_sample
# write.table(allpop_sam, file="Z:/allpops_sample.txt",row.names=F,col.names=F,quote=F,sep="\t")
# 
# 
# #make spoofed.fam for plink2
# #requires 0 sampleid 0 0 0 0
# 
# spoofed <- data.frame(id1=rep(0,nrow(allpop_sam)), id2=rep(0,nrow(allpop_sam)), id3=rep(0,nrow(allpop_sam)),
#                       id2=rep(0,nrow(allpop_sam)))
# spoofed <- cbind(allpop_sam, spoofed)
# write.table(spoofed, file="Z:/data/spoofed.fam", row.names=F,col.names=F,quote=F,sep="\t")
