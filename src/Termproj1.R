getwd()

setwd("~/Desktop/SNU/2021-1/Bioinfo1/dat")

#Data loading
df<-read.table("read-counts.txt")
head(df)

colnames(df)<-df[1,]
df<-df[-1,]
head(df)

#1. Original scenario
dfs<-df

CLIP<-as.numeric(dfs$`CLIP-35L33G.bam`)
RNA.control<-as.numeric(dfs$`RNA-control.bam`)
chip.enrichment<-CLIP/RNA.control

RPF<-as.numeric(dfs$`RPF-siLin28a.bam`)
RNA<-as.numeric(dfs$`RNA-siLin28a.bam`)
rden.change<-RPF/RNA

#Data preprocessing
for.plot <- data.frame(log2(chip.enrichment),log2(rden.change))

for.plot<-na.omit(for.plot) #Removing Nan index
for.plot<-for.plot[is.finite(rowSums(for.plot)),] #Removing inf index
dim(for.plot)

#Scatter plotting
plot(for.plot[,1],for.plot[,2]) #The shape of plot is not similar to the example...

#Checking Spearman Correlation
##Spearman correlation was used to be less affected by outliers.
cor(for.plot[,1],for.plot[,2],method="spearman") # 0.13 Too low :(

#Check the distributed length
len <- as.numeric(df$Length)
mean(len)
median(len)
boxplot(len) #I think data need to be Filtered or Normalized!!

#2-1. Filtered scenario
#Filtering with arbitrary threshold.
library(dplyr)

thresh <- 5000
df2<-filter(df, Length > thresh)

dfs<-df2

CLIP<-as.numeric(dfs$`CLIP-35L33G.bam`)
RNA.control<-as.numeric(dfs$`RNA-control.bam`)
chip.enrichment<-CLIP/RNA.control

RPF<-as.numeric(dfs$`RPF-siLin28a.bam`)
RNA<-as.numeric(dfs$`RNA-siLin28a.bam`)
rden.change<-RPF/RNA

#Data preprocessing
for.plot <- data.frame(log2(chip.enrichment),log2(rden.change))

for.plot<-na.omit(for.plot) #Removing Nan index
for.plot<-for.plot[is.finite(rowSums(for.plot)),] #Removing inf index
dim(for.plot)

#Scatter plotting
plot(for.plot[,1],for.plot[,2])

#Checking Spearman Correlation
cor(for.plot[,1],for.plot[,2],method="spearman") # 0.157 (>0.13) 

#2-2. Filtered scenario
#Checking trend with thresholds.
seq.Thresh <- seq(from=0,to=20000, by=100)

my.cor <- function(df,thresh,method="spearman"){
  df2<-filter(df, Length > thresh)
  
  dfs<-df2
  
  CLIP<-as.numeric(dfs$`CLIP-35L33G.bam`)
  RNA.control<-as.numeric(dfs$`RNA-control.bam`)
  chip.enrichment<-CLIP/RNA.control
  
  RPF<-as.numeric(dfs$`RPF-siLin28a.bam`)
  RNA<-as.numeric(dfs$`RNA-siLin28a.bam`)
  rden.change<-RPF/RNA
  
  #Data preprocessing
  for.plot <- data.frame(log2(chip.enrichment),log2(rden.change))
  
  for.plot<-na.omit(for.plot) #Removing Nan index
  for.plot<-for.plot[is.finite(rowSums(for.plot)),] #Removing inf index
  
  return(cor(for.plot[,1],for.plot[,2],method=method))
}

cor.vec <- c()
for(trsh in seq.Thresh){
  cor.vec<-c(cor.vec,my.cor(df,trsh))
}
cor.vec

plot(seq.Thresh,cor.vec)

#Trend for the Spearman correlation
pearson.cor.vec <- c()
for(trsh in seq.Thresh){
  pearson.cor.vec<-c(pearson.cor.vec,my.cor(df,trsh,method="pearson"))
}
pearson.cor.vec

plot(seq.Thresh,pearson.cor.vec)

#Just for comparison between Pearson and Spearman
summary(pearson.cor.vec)
summary(cor.vec)

#If there is no threshold from any recognized value for length,
#I would like to use threshold near 5000.

#Result plot with type label.
##Used Threshold: 5000
thresh <- 5000
df2<-filter(df, Length > thresh)

ref<-read.csv(url('https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt'), sep='\t')
colnames(ref)[1]<-"Geneid"

#Preprocessing
df2$Geneid<-substr(df2$Geneid,1,18)

dfs<-merge(df2,ref,by="Geneid")
colnames(dfs)

CLIP<-as.numeric(dfs$`CLIP-35L33G.bam`)
RNA.control<-as.numeric(dfs$`RNA-control.bam`)
chip.enrichment<-CLIP/RNA.control

RPF<-as.numeric(dfs$`RPF-siLin28a.bam`)
RNA<-as.numeric(dfs$`RNA-siLin28a.bam`)
rden.change<-RPF/RNA

#Data preprocessing
for.plot <- data.frame(log2(chip.enrichment),log2(rden.change),dfs$type)

for.plot<-na.omit(for.plot) #Removing Nan index
for.plot<-for.plot[is.finite(rowSums(for.plot[,1:2])),] #Removing inf index
dim(for.plot)
colnames(for.plot)<-c("log2_Chip.Enrichment","log2_rden.change","Type")

library(ggplot2)

png(filename="../rst/myplot.png",width=1400,height=1200)
ggplot(data=for.plot, mapping=aes(x=log2_Chip.Enrichment,y=log2_rden.change)) +
  geom_point(aes(color=Type))
dev.off()

cor(for.plot$log2_Chip.Enrichment,for.plot$log2_rden.change,method = "spearman") # 0.28
cor(for.plot$log2_Chip.Enrichment,for.plot$log2_rden.change,method = "pearson") # 0.34 Wow

#The correlations from Original filtered set and from genes with type informations set are quite different.
#I think the set which has type information cannot represent the whole data set well.

