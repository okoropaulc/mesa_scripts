model <- c("EN", "EN", "EN", "RF", "RF", "RF", "SVR", "SVR", "SVR", "KNN", "KNN", "KNN")
threshold <- c("all","p>0", "p>0.1")


afa_count <- c(3268, 2351, 1545, 2673, 1858, 1167, 2492, 1663, 961, 2333, 1500, 824)
afa <- data.frame(population=rep("AFA", 12), Model=model, Threshold=rep(threshold, 4), Count=afa_count)

cau_count <- c(4266, 2659, 1447, 2915, 1845, 1027, 3188, 1994, 1058, 2457, 1459, 748)
cau <- data.frame(population=rep("CAU", 12), Model=model, Threshold=rep(threshold, 4), Count=cau_count)

his_count <- c(3911, 2527, 1505, 2804, 1912, 1145, 2882, 1905, 1051, 2389, 1548, 786)
his <- data.frame(population=rep("HIS", 12), Model=model, Threshold=rep(threshold, 4), Count=his_count)

afhi_count <- c(4975, 3483, 2150, 3375, 2335, 1450, 3345, 2280, 1364, 2622, 1745, 968)
afhi <- data.frame(population=rep("AFHI", 12), Model=model, Threshold=rep(threshold, 4), Count=afhi_count)

all_count <- c(5231, 3638, 2218, 3429, 2365, 1418, 3553, 2399, 1359, 2457, 1617, 865)
all <- data.frame(population=rep("ALL", 12), Model=model, Threshold=rep(threshold, 4), Count=all_count)

table2 <- rbind(afa, cau, his, afhi, all)

library(ggplot2)
ggplot(data=table2, aes(x=population, y=Count, fill=Model, shape=Threshold)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal(30) + xlab("Population") + ylab("Genes")

ggplot(data=table2, aes(x=population, y=Count, fill=Model, color=Threshold)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal(30)  + xlab("Population") + ylab("Genes")

ggplot(data=table2, aes(x=population, y=Count, fill=Model, color=Threshold)) +
  geom_bar(stat="identity", position=position_dodge2())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()  + xlab("Population") + ylab("Genes")

library(viridis)

ggplot(data=table2, aes(x=population, y=Count, color=Model, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") +
  theme_minimal(30)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T)

ggplot(data=table2, aes(x=population, y=Count, colour=Model, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", aesthetics = "colour") +
  theme_minimal(30)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T, option="A")

ggplot(data=table2, aes(x=population, y=Count, colour=Threshold, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", aesthetics = "fill") +
  theme_minimal(30)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T, option="A")


ggplot(data=table2, aes(x=population, y=Count, colour=Model, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", aesthetics = "colour") +
  theme_minimal(30)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T, option="A")+
  geom_text(aes(label=Count), vjust=2, color="white", size=4)


#face wrap

ggplot(data=table2, aes(x=population, y=Count, colour=Model, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", aesthetics = "colour") +
  theme_minimal(30)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T, option="A")+
  facet_wrap(~Threshold)


ggplot(data=table2, aes(x=population, y=Count, colour=Threshold, fill=Threshold)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", aesthetics = "colour") +
  theme_minimal(30)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T, option="D")+
  facet_wrap(~Model)
#width=1300 height=1100


ggplot(data=table2, aes(x=population, y=Count, fill=Threshold)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", aesthetics = "colour") +
  theme_minimal(30)  + xlab("Population") + ylab("Genes") + scale_fill_viridis(discrete=T, option="A")+
  facet_wrap(~Model) + scale_fill_discrete(name =expression(rho), labels = c("all", ">0", ">0.1"))
#width=1300 height=1100