#METS predicting Expression in MESA
library(stringr)
#CAU
mets2cau <- read.table("/home/paul/METS_model/METS_to_MESA_CAU_unfiltered_models_predicted_expression.txt", header = TRUE)
mets2cau$FID <- NULL 
rownames(mets2cau) <-mets2cau$IID
mets2cau$IID <- NULL
for (i in 1:length(mets2cau)){
  colnames(mets2cau)[i] <- gsub('\\.[0-9]+','',names(mets2cau[i]))
}
probeid_m2cau <- colnames(mets2cau)
sampleid_m2cau <- rownames(mets2cau)

caugex <- read.table(file = "/home/paul/METS_model/CAU_PF10.txt", header = T)
for (i in 2:length(caugex)){
  colnames(caugex)[i] <- str_sub(names(caugex[i]),2,6)
}
caugex$PROBE_ID <- as.character(caugex$PROBE_ID)
for (i in 1:length(caugex$PROBE_ID)){
  caugex$PROBE_ID[i] <- gsub('\\.[0-9]+','',caugex$PROBE_ID[i])
}
caugex <- caugex[which(caugex$PROBE_ID %in% names(mets2cau)),]
caugex <- caugex[order(caugex$PROBE_ID),] #order the dataframe by probe ID
rownames(caugex) <- c(1:length(caugex$PROBE_ID))

mets2cau <- as.data.frame(t(mets2cau))
mets2cau <- cbind(PROBE_ID = probeid_m2cau, mets2cau)
mets2cau <- mets2cau[which(mets2cau$PROBE_ID %in% caugex$PROBE_ID),]
mets2cau <- mets2cau[order(mets2cau$PROBE_ID),]
rownames(mets2cau) <- c(1:length(mets2cau$PROBE_ID)) #rename the rows to numbers

mets2cau <- mets2cau %>% dplyr::select(names(caugex))

m2cau <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(caugex$PROBE_ID)){
  m2cau[i,1] <- caugex$PROBE_ID[i]
  m2cau[i,2] <- cor(as.numeric(caugex[i,2:length(caugex)]), as.numeric(mets2cau[i,2:length(mets2cau)]), method = "spearman")
}
colnames(m2cau) <- c("gene", "spearman")
library(tidyr)
m2cau <- m2cau %>% drop_na() #drop rows containing NA
write.table(m2cau, file = "/home/paul/METS_model/METS_2_CAU_spearman.txt", row.names = FALSE) #write out results



#HIS
mets2his <- read.table("/home/paul/METS_model/METS_to_MESA_HIS_unfiltered_models_predicted_expression.txt", header = TRUE)
mets2his$FID <- NULL 
rownames(mets2his) <-mets2his$IID
mets2his$IID <- NULL
for (i in 1:length(mets2his)){
  colnames(mets2his)[i] <- gsub('\\.[0-9]+','',names(mets2his[i]))
}
probeid_m2his <- colnames(mets2his)
sampleid_m2his <- rownames(mets2his)

hisgex <- read.table(file = "/home/paul/METS_model/HIS_PF10.txt", header = T)
for (i in 2:length(hisgex)){
  colnames(hisgex)[i] <- str_sub(names(hisgex[i]),2,6)
}
hisgex$PROBE_ID <- as.character(hisgex$PROBE_ID)
for (i in 1:length(hisgex$PROBE_ID)){
  hisgex$PROBE_ID[i] <- gsub('\\.[0-9]+','',hisgex$PROBE_ID[i])
}
hisgex <- hisgex[which(hisgex$PROBE_ID %in% names(mets2his)),]
hisgex <- hisgex[order(hisgex$PROBE_ID),] #order the dataframe by probe ID
rownames(hisgex) <- c(1:length(hisgex$PROBE_ID))

mets2his <- as.data.frame(t(mets2his))
mets2his <- cbind(PROBE_ID = probeid_m2his, mets2his)
mets2his <- mets2his[which(mets2his$PROBE_ID %in% hisgex$PROBE_ID),]
mets2his <- mets2his[order(mets2his$PROBE_ID),]
rownames(mets2his) <- c(1:length(mets2his$PROBE_ID)) #rename the rows to numbers

mets2his <- mets2his %>% dplyr::select(names(hisgex))

m2his <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(hisgex$PROBE_ID)){
  m2his[i,1] <- hisgex$PROBE_ID[i]
  m2his[i,2] <- cor(as.numeric(hisgex[i,2:length(hisgex)]), as.numeric(mets2his[i,2:length(mets2his)]), method = "spearman")
}
colnames(m2his) <- c("gene", "spearman")
library(tidyr)
m2his <- m2his %>% drop_na() #drop rows containing NA
write.table(m2his, file = "/home/paul/METS_model/METS_2_HIS_spearman.txt", row.names = FALSE) #write out results


#AFHI
mets2afhi <- read.table("/home/paul/METS_model/METS_to_MESA_AFHI_unfiltered_models_predicted_expression.txt", header = TRUE)
mets2afhi$FID <- NULL 
rownames(mets2afhi) <-mets2afhi$IID
mets2afhi$IID <- NULL
for (i in 1:length(mets2afhi)){
  colnames(mets2afhi)[i] <- gsub('\\.[0-9]+','',names(mets2afhi[i]))
}
probeid_m2afhi <- colnames(mets2afhi)
sampleid_m2afhi <- rownames(mets2afhi)

afhigex <- read.table(file = "/home/paul/METS_model/AFHI_PF10.txt", header = T)
for (i in 2:length(afhigex)){
  colnames(afhigex)[i] <- str_sub(names(afhigex[i]),2,6)
}


afhigex$id <- as.character(afhigex$id)
for (i in 1:length(afhigex$id)){
  afhigex$id[i] <- gsub('\\.[0-9]+','',afhigex$id[i])
}
afhigex <- afhigex[which(afhigex$id %in% names(mets2afhi)),]
afhigex <- afhigex[order(afhigex$id),] #order the dataframe by probe ID
rownames(afhigex) <- c(1:length(afhigex$id))

mets2afhi <- as.data.frame(t(mets2afhi))
mets2afhi <- cbind(PROBE_ID = probeid_m2afhi, mets2afhi)
mets2afhi <- mets2afhi[which(mets2afhi$PROBE_ID %in% afhigex$id),]
mets2afhi <- mets2afhi[order(mets2afhi$PROBE_ID),]
rownames(mets2afhi) <- c(1:length(mets2afhi$PROBE_ID))#rename the rows to numbers
colnames(mets2afhi)[1] <- "id"
mets2afhi <- mets2afhi %>% dplyr::select(names(afhigex))

m2afhi <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(afhigex$id)){
  m2afhi[i,1] <- afhigex$id[i]
  m2afhi[i,2] <- cor(as.numeric(afhigex[i,2:length(afhigex)]), as.numeric(mets2afhi[i,2:length(mets2afhi)]), method = "spearman")
}
colnames(m2afhi) <- c("gene", "spearman")
library(tidyr)
m2afhi <- m2afhi %>% drop_na() #drop rows containing NA
write.table(m2afhi, file = "/home/paul/METS_model/METS_2_AFHI_spearman.txt", row.names = FALSE) #write out results


#ALL
mets2all <- read.table("/home/paul/METS_model/METS_to_MESA_ALL_unfiltered_models_predicted_expression.txt", header = TRUE)
mets2all$FID <- NULL 
rownames(mets2all) <-mets2all$IID
mets2all$IID <- NULL
for (i in 1:length(mets2all)){
  colnames(mets2all)[i] <- gsub('\\.[0-9]+','',names(mets2all[i]))
}
probeid_m2all <- colnames(mets2all)
sampleid_m2all <- rownames(mets2all)

allgex <- read.table(file = "/home/paul/METS_model/ALL_PF10.txt", header = T)
for (i in 2:length(allgex)){
  colnames(allgex)[i] <- str_sub(names(allgex[i]),2,6)
}


allgex$id <- as.character(allgex$id)
for (i in 1:length(allgex$id)){
  allgex$id[i] <- gsub('\\.[0-9]+','',allgex$id[i])
}
allgex <- allgex[which(allgex$id %in% names(mets2all)),]
allgex <- allgex[order(allgex$id),] #order the dataframe by probe ID
rownames(allgex) <- c(1:length(allgex$id))

mets2all <- as.data.frame(t(mets2all))
mets2all <- cbind(PROBE_ID = probeid_m2all, mets2all)
mets2all <- mets2all[which(mets2all$PROBE_ID %in% allgex$id),]
mets2all <- mets2all[order(mets2all$PROBE_ID),]
rownames(mets2all) <- c(1:length(mets2all$PROBE_ID))#rename the rows to numbers
colnames(mets2all)[1] <- "id"
mets2all <- mets2all %>% dplyr::select(names(allgex))

m2all <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(allgex$id)){
  m2all[i,1] <- allgex$id[i]
  m2all[i,2] <- cor(as.numeric(allgex[i,2:length(allgex)]), as.numeric(mets2all[i,2:length(mets2all)]), method = "spearman")#calculate spearman correlation between measured and predicted gene expression
}
colnames(m2all) <- c("gene", "spearman")
library(tidyr)
m2all <- m2all %>% drop_na() #drop rows containing NA
write.table(m2all, file = "/home/paul/METS_model/METS_2_ALL_spearman.txt", row.names = FALSE) #write out results



#AFA
mets2afa <- read.table("/home/paul/METS_model/METS_to_MESA_AFA_unfiltered_models_predicted_expression.txt", header = TRUE)
mets2afa$FID <- NULL 
rownames(mets2afa) <-mets2afa$IID
mets2afa$IID <- NULL
for (i in 1:length(mets2afa)){
  colnames(mets2afa)[i] <- gsub('\\.[0-9]+','',names(mets2afa[i]))
}
probeid_m2afa <- colnames(mets2afa)
sampleid_m2afa <- rownames(mets2afa)

afagex <- read.table(file = "/home/paul/METS_model/AFA_PF10.txt", header = T)
for (i in 2:length(afagex)){
  colnames(afagex)[i] <- str_sub(names(afagex[i]),2,6)
}
afagex$PROBE_ID <- as.character(afagex$PROBE_ID)
for (i in 1:length(afagex$PROBE_ID)){
  afagex$PROBE_ID[i] <- gsub('\\.[0-9]+','',afagex$PROBE_ID[i])
}
afagex <- afagex[which(afagex$PROBE_ID %in% names(mets2afa)),]
afagex <- afagex[order(afagex$PROBE_ID),] #order the dataframe by probe ID
rownames(afagex) <- c(1:length(afagex$PROBE_ID))

mets2afa <- as.data.frame(t(mets2afa))
mets2afa <- cbind(PROBE_ID = probeid_m2afa, mets2afa)
mets2afa <- mets2afa[which(mets2afa$PROBE_ID %in% afagex$PROBE_ID),]
mets2afa <- mets2afa[order(mets2afa$PROBE_ID),]
rownames(mets2afa) <- c(1:length(mets2afa$PROBE_ID)) #rename the rows to numbers

mets2afa <- mets2afa %>% dplyr::select(names(afagex))

m2afa <- data.frame(gene = NULL, spearman = NULL) #create empty dataframe to store results
for (i in 1:length(afagex$PROBE_ID)){
  m2afa[i,1] <- afagex$PROBE_ID[i]
  m2afa[i,2] <- cor(as.numeric(afagex[i,2:length(afagex)]), as.numeric(mets2afa[i,2:length(mets2afa)]), method = "spearman")#calculate spearman correlation between measured and predicted gene expression
}
colnames(m2afa) <- c("gene", "spearman")
library(tidyr)
m2afa <- m2afa %>% drop_na() #drop rows containing NA
write.table(m2afa, file = "/home/paul/METS_model/METS_2_AFA_spearman.txt", row.names = FALSE) #write out results


#plotting the results in a violin plot
#my expression data
m2CAU <- read.table(file = "Z:/METS_model/METS_2_CAU_spearman.txt", header = T)
m2c <- subset(m2CAU, m2CAU$spearman > 0.1)
ALL <- read.table(file = "Z:/METS_model/METS_2_ALL_spearman.txt", header = T)
m2HIS <- read.table(file = "Z:/METS_model/METS_2_HIS_spearman.txt", header = T)
m2h <- subset(m2HIS, m2HIS$spearman > 0.1)
m2AFA <- read.table(file = "Z:/METS_model/METS_2_AFA_spearman.txt", header = T)
m2a <- subset(m2AFA, m2AFA$spearman > 0.1)
AFHI <- read.table(file = "/home/paul/METS_model/METS_2_AFHI_spearman.txt", header = T)

library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

# colnames(ALL)<-c("gene.ALL","spearman.ALL")
# colnames(AFHI)<-c("gene.AFHI","spearman.AFHI")
# colnames(AFA)<-c("gene.AFA","spearman.AFA")
# colnames(CAU)<-c("gene.CAU","spearman.CAU")
# colnames(HIS)<-c("gene.HIS","spearman.HIS")

ALL$pop<-"ALL"
ALL$median<-round(median(ALL$spearman),3)
AFHI$pop<-"AFHI"
AFHI$median<-round(median(AFHI$spearman),3)
AFA$pop<-"AFA"
AFA$median<-round(median(AFA$spearman),3)
CAU$pop<-"CAU"
CAU$median<-round(median(CAU$spearman),3)
HIS$pop<-"HIS"
HIS$median<-round(median(HIS$spearman),3)
total<-rbind.data.frame(ALL,AFHI,AFA,CAU,HIS)
total[is.na(total)]<-0
pV<-ggplot(data=total,aes(y=spearman,x=as.factor(median)))
pV + geom_violin(draw_quantiles = T, aes(fill=pop)) + 
  geom_boxplot(width = 0.4) +
  geom_hline(yintercept=0) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("MESA spearman correlation using METS models chr 1, 9 - 22")

#ryan expression data

CAU <- read.table(file = "Z:/METS_model/METS_2_CAU_spearman2.txt", header = T)
ALL <- read.table(file = "Z:/METS_model/METS_2_ALL_spearman2.txt", header = T)
HIS <- read.table(file = "Z:/METS_model/METS_2_HIS_spearman2.txt", header = T)
AFA <- read.table(file = "Z:/METS_model/METS_2_AFA_spearman2.txt", header = T)
AFHI <- read.table(file = "Z:/METS_model/METS_2_AFHI_spearman2.txt", header = T)

library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

# colnames(ALL)<-c("gene.ALL","spearman.ALL")
# colnames(AFHI)<-c("gene.AFHI","spearman.AFHI")
# colnames(AFA)<-c("gene.AFA","spearman.AFA")
# colnames(CAU)<-c("gene.CAU","spearman.CAU")
# colnames(HIS)<-c("gene.HIS","spearman.HIS")

ALL$pop<-"ALL"
ALL$median<-round(median(ALL$spearman),3)
AFHI$pop<-"AFHI"
AFHI$median<-round(median(AFHI$spearman),3)
AFA$pop<-"AFA"
AFA$median<-round(median(AFA$spearman),3)
CAU$pop<-"CAU"
CAU$median<-round(median(CAU$spearman),3)
HIS$pop<-"HIS"
HIS$median<-round(median(HIS$spearman),3)
total<-rbind.data.frame(ALL,AFHI,AFA,CAU,HIS)
total[is.na(total)]<-0
pV<-ggplot(data=total,aes(y=spearman,x=as.factor(median)))
pV + geom_violin(draw_quantiles = T, aes(fill=pop)) + 
  geom_boxplot(width = 0.4) +
  geom_hline(yintercept=0) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("MESA spearman correlation using METS models chr 1, 9 - 22")

ptm <- proc.time()
print("yes")
ptm2 <- proc.time() - ptm

fileConn<-file("/home/paul/mesa_models/svr/try.txt")
for (i in 1:4){
  ptm <- proc.time()
  pV + geom_violin(draw_quantiles = T, aes(fill=pop)) + 
    geom_boxplot(width = 0.4) +
    geom_hline(yintercept=0) +
    scale_color_viridis(discrete=TRUE) +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("MESA spearman correlation using METS models chr 1, 9 - 22")
  ptm2 <- proc.time() - ptm
  ptm2 <- as.matrix(ptm2)
  write(ptm2[3], fileConn, append = TRUE, sep = "\n")
}
close(fileConn)

