#Make Figure 6
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
  ifelse(is.na(en_rf[i,2]) | is.na(en_rf[i,3]), en_rf$type[i]<-"unique", en_rf$type[i]<-"common")
}

#now replace all NAs with 0
en_rf[is.na(en_rf)] <- 0

#plot them
ggplot(en_rf, aes(x=en, y=rf, col=type)) + geom_point(lwd=3)  +
  xlab("Elastic Net") + ylab("Random Forest") +
  theme_classic(30) + geom_abline(intercept=0, slope=1, color="blue") +
  scale_color_manual(values = c("black","red"))#  + ggtitle("EN vs RF on HDL")
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