#
## Authors= c("M. Liu","J Grayson", "S Simpson) 




## Call this function in the console to import necessary packages
## Necessary if running on new computer with fresh R installation

ImportPackages<- function(){
  source("https://bioconductor.org/biocLite.R")
  biocLite("flowCore")
  biocLite("flowClust")
  biocLite("flowTrans")
  
  install.packages("ggplot2")
  install.packages("tidyr")
  install.packages("plyr")
  install.packages("RColorBrewer")
}

## Clean up workspace
rm(list = ls())

## Bring in some libraries to read in FCS files
library(flowCore)
library(flowClust)
library(flowTrans)

# set working directory
setwd()

## Prompt the user to provide a desired number of clusters to break the data into
##
## Code needs to be sourced for this to work
##
## k<-as.numeric(readline("How many clusters?"))

## Define the number of clusters to be used in the analysis
k<-6

# list all files from the current directory.Use the pattern argument to define a common pattern  for import files with regex. Here: .fcs
list.files(pattern=".fcs$") 

## Read in tab delimited file containing mouse labels (column 1) and
## treatment labels (column 2) WITH COLUMN HEADERS
##
## This could be expanded to bring in a list of markers and their associated fluorescence
## parameters
Label <- data.frame(read.table("Mouse label table.txt", header=T))
Label

## create a list from these files
list.filenames<-list.files(pattern=".fcs$")
print(list.filenames)

## create an empty list that will serve as a container to receive the incoming files
list.data<-list()

# create a loop to read in your data
for (i in 1:length(list.filenames))
{
    list.data[[i]]<-read.FCS(list.filenames[[i]],transformation = FALSE)
    list.data[[i]]<-as.data.frame(exprs(list.data[[i]]))
    ## take out FSC,SSC,TIME
    list.data[[i]] <- list.data[[i]][,c(-1,-2,-3)]
    ## remove saturated values (defined manually)
    RM<-(list.data[[i]][,1]< -100)| (list.data[[i]][,1] >40000)|
        (list.data[[i]][,2]< -100)| (list.data[[i]][,2]> 20000)|
        (list.data[[i]][,3]< -100)| (list.data[[i]][,3]> 20000)|
        (list.data[[i]][,4]< -100)| (list.data[[i]][,4]> 10000)|
        (list.data[[i]][,5]< -100)| (list.data[[i]][,5]> 40000)|
        (list.data[[i]][,6]< -100)| (list.data[[i]][,6]> 30000)|
        (list.data[[i]][,7]< -100)| (list.data[[i]][,7]> 100000)
    list.data[[i]]<-list.data[[i]][!RM,]
    ##remove variables not included in the model
    list.data[[i]] <- list.data[[i]][,c(-8)]
    ##add labels to dataset
    list.data[[i]] <- cbind(list.data[[i]], Mouse = rep(Label$Mouse[i],nrow(list.data[[i]])), Treatment = rep(Label$Treatment[i],nrow(list.data[[i]])))
}

##convert the list to a dataframe
s5_with_label <- NULL
for (i in 1:length(list.filenames))
{
    s5_with_label<- rbind(s5_with_label,list.data[[i]])
}

##name the variables
colnames(s5_with_label) <- c("","mouse","treatment") 

##save the data
save(s5_with_label,file="")




##Install packages
library(flowCore)
library(flowTrans) 
library(ggplot2) 
library(tidyr)
library(plyr)
library(RColorBrewer)


##read the written Rda file
load("") 
##take out mouse label and treatment type for data transformation
ALL<-s5_with_label[,c(-8,-9)] 

##ArcSinh transformation using flowTrans package
ALL<-as.matrix(ALL)

## Convert matrix to flow frame and perform ArcSin transformation
ALL<-flowFrame(ALL) 
ALL_trans<-flowTrans(dat=ALL,fun="mclMultivArcSinh",colnames(ALL), n2f=FALSE,parameters.only = FALSE)

## Extract the transformed data
ALL_after_transformation<-exprs(ALL_trans$result)  

## Scale the transformed data
ALL<-as.data.frame(scale(ALL_after_transformation)) 

## Determine optimal number of clusters by WSS plot
set.seed(2L)
wss <- (nrow(ALL)-1)*sum(apply(ALL,2,var))

#plot k by wss to find the elbow of the plot
for (i in 2:25) wss[i] <- sum(kmeans(ALL,algorithm = "Lloyd", centers=i)$withinss)
plot(1:25, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

## Prompt the user to provide a desired number of clusters to break the data into
##
## Code needs to be sourced for this to work
##
## k<-as.numeric(readline("How many clusters?"))

## Define the number of clusters to be used in the analysis
k<-6

##K-means clustering
set.seed(1L)  
ALL_cluster <- kmeans(ALL,k,algorithm = "Lloyd", nstart = 1000)

#Add cluster number to the dataset
ALL_after_cluster <- cbind(ALL,cluster=ALL_cluster$cluster) 

##Visualize the characteristics of clusters
##stacked density plot to show the characteristic of each cluster



##calculate the density of fluorescent molecules depending on the cluster
Density=NULL
for (i in 1:7)
{
    #take out cluster and combine with each fluorscence in each round
    my.data <- as.data.frame(cbind( V1 = ALL_after_cluster[,i], cluster=ALL_after_cluster[,ncol(ALL_after_cluster)])) 
    
    library(plyr)
    ##calculate density of chosen molecules depending on cluster
    res <- dlply(my.data, .(cluster), function(x) density(x$V1)) 
    ##take out the density and value from the result list
    dd <- ldply(res, function(z){data.frame(Values = z[["x"]], 
                                            V1_density = z[["y"]],
                                            V1_count = z[["y"]]*z[["n"]])}) 
    #adjust the offset
    dd$offest <- -0.4*dd$cluster 
    ##add the offset to the density to avoid verical overlap
    dd$V1_density_offest=dd$V1_density+dd$offest 
    Fluor <- colnames(ALL_after_cluster)[i] 
    
    ##generate a dataset combining offset, density and original value
    DD <- cbind(dd$cluster,dd$Values,dd$V1_density_offest,dd$offest) 
    colnames(DD) <- c("cluster",paste(Fluor,"Value"),paste(Fluor,"density_offset"),paste(Fluor,"offset"))
    #combine results of different molecules
    Density <- cbind(Density,DD) 
}


##separate the dataset by molecules and add label
CD3_density <- as.data.frame(Density[,1:4])
colnames(CD3_density)<-c("cluster","FI","density","offset")

CD69_density <- as.data.frame(Density[,5:8])
colnames(CD69_density)<-c("cluster","FI","density","offset")

CD62L_density <- as.data.frame(Density[,9:12])
colnames(CD62L_density)<-c("cluster","FI","density","offset")

CD90.2_density <- as.data.frame(Density[,13:16])
colnames(CD90.2_density)<-c("cluster","FI","density","offset")

CD25_density <- as.data.frame(Density[,17:20])
colnames(CD25_density)<-c("cluster","FI","density","offset")

CD4_density <- as.data.frame(Density[,21:24])
colnames(CD4_density)<-c("cluster","FI","density","offset")

CD44_density <- as.data.frame(Density[,25:28])
colnames(CD44_density)<-c("cluster","FI","density","offset")

#Converts name of an object to a string and returns that string
deparser <- function(y){
  
  deparsedStr<-deparse(substitute(y))
  return(deparsedStr)
  
}

# Exports plot as PNG file
exportPNG<- function(imagePlot){
  
  imagePlotName<-paste(deparser(imagePlot),".png")
  png(imagePlotName)
  plot(imagePlot)
  dev.off()
  
}

# Function to generate density plots for each cluster for a given parameter
# Exports graph as a PNG to the working directory
markerHistograms<-function(dFrame){
  
  dPlot<-ggplot(dFrame, aes(x=FI, y=density,group=cluster)) + geom_line(color="white") + 
  geom_ribbon(aes(FI, ymin=offset,ymax=density),alpha=0.5) + scale_y_continuous(breaks=NULL) + 
  ggtitle("") + xlab("") + ylab("") + theme(panel.background = element_rect(fill="white",colour = "black", size = 0.75), 
  axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(),axis.ticks=element_blank()) + 
  theme_bw() + ggtitle(deparse(substitute(dFrame))) + theme(plot.title = element_text(hjust = 0.5))
  
  imagePlotName<-paste(deparse(substitute(dFrame)),".png")
  png(imagePlotName)
  plot(dPlot)
  dev.off()
  
  return(dPlot)
}

markerHistograms(CD3_density)
markerHistograms(CD69_density)
markerHistograms(CD62L_density)
markerHistograms(CD90.2_density)
markerHistograms(CD25_density)
markerHistograms(CD4_density)
markerHistograms(CD44_density)

##add mouse label to the dataset
ALL_after_cluster<- cbind(ALL_after_cluster, mouse=as.matrix(s5_with_label)[,ncol(s5_with_label)-1])

##create a function to calculate percentage of each cluster
r=0
cluster_ratio<-function(x){
    for (i in 1:k){
        count<-length(x[x==i])
        ratio<-round(count/length(x),4)*100
        r[i]=ratio
    }
    return(r)
}

##apply the cluster_ratio function to calculate the distribution of cell to each cluster
result_of_all<-as.matrix(aggregate(ALL_after_cluster$cluster, by=list(ALL_after_cluster$mouse), FUN=cluster_ratio))

### Rename the dataset

## Generate a vector containing Cluster labels, depending on the the number of clusters
## used in the analysis, e.g. "Cluster 1, Cluster 2, Cluster 3..."

ClustNum<-c(1:k)
ClustLabels<-NULL

for(i in 1:k){
  
ClustLabels<-c(ClustLabels, paste("Cluster",as.character(ClustNum[i])))

}

## Apply cluster labels to results data frame
colnames(result_of_all)<-c("mouse", ClustLabels)
print(result_of_all)

## create treatment label
Treatment<-as.character(Label$Treatment)

#add the label to result_of_all
result_of_all<-as.data.frame(cbind(result_of_all,Treatment))
print(result_of_all)

## write as CSV
## write.csv(result_of_all, file = "CD4 Cluster.csv", append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE,
##          col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

##remove treatment label
result_of_all_1 <- result_of_all[,c(-1,-(k+2))]

##create a convert_to_numeric function
numeric_conversion<- function(x){
    for (i in 1:ncol(x)){
        x[,i] = as.numeric(as.character(x[,i]))
    }
    return(x)
}

##convert factor in result_of_all_1 to numeric variables
result_of_all <- cbind(numeric_conversion(result_of_all_1),treatment = result_of_all[,k+2])

##Aggregate the result by infections for analysis
##gather the data using tidyr package
ALL_gather <- gather(result_of_all,"cluster","percentage_of_cluster", -treatment) 

##calculate the mean percentage of cells in each cluster
mean_for_rose<-as.data.frame(aggregate(ALL_gather$percentage_of_cluster, by=list(ALL_gather$treatment,ALL_gather$cluster), FUN=mean))

##calculate the standard deviation
sd_for_rose<-as.data.frame(aggregate(ALL_gather$percentage_of_cluster, by=list(ALL_gather$treatment,ALL_gather$cluster), FUN=sd))

##generate a dataset that have mean and sd
ALL_for_rose <-cbind(mean_for_rose,sd_for_rose[,3])
colnames(ALL_for_rose) <- c("treatment","cluster","percentage_of_cluster","sd")

##generate a barplot showing the distribution of CD4 cells to the 12 clusters
##plot all treatments
ALL_bar <- ggplot(ALL_for_rose, aes(x=cluster, y=percentage_of_cluster,group=cluster, fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd),colour="black") +facet_grid(.~treatment)+ scale_fill_brewer(palette = "Set3") + scale_color_brewer()

png("Bar graphs, all groups.png")
plot(ALL_bar)
dev.off()

#plot each treatment individually


# D18N_bar<- ggplot(ALL_for_rose[ALL_for_rose$treatment=="D18N",], aes(x=cluster, y=percentage_of_cluster, group= cluster,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer(palette = "Set3") + scale_color_brewer()
# 
# WT_bar <- ggplot(ALL_for_rose[ALL_for_rose$treatment=="WT",], aes(x=cluster, y=percentage_of_cluster,group=cluster, fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd),colour = "black")+ scale_fill_brewer(palette = "Set3") + scale_color_brewer()
# 
# DNDN_bar<- ggplot(ALL_for_rose[ALL_for_rose$treatment=="D18N>D18N",], aes(x=cluster, y=percentage_of_cluster, group= cluster,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer(palette = "Set3") + scale_color_brewer()
# 
# DNWT_bar<- ggplot(ALL_for_rose[ALL_for_rose$treatment=="D18N>WT",], aes(x=cluster, y=percentage_of_cluster, group= cluster,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer(palette = "Set3") + scale_color_brewer()
# 
# WTWT_bar<- ggplot(ALL_for_rose[ALL_for_rose$treatment=="WT>WT",], aes(x=cluster, y=percentage_of_cluster, group= cluster,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer(palette = "Set3") + scale_color_brewer()
# 
# WTDN_bar<- ggplot(ALL_for_rose[ALL_for_rose$treatment=="WT>DN",], aes(x=cluster, y=percentage_of_cluster, group= cluster,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer(palette = "Set3") + scale_color_brewer()

# WT_bar
# D18N_bar
# WTWT_bar
# WTDN_bar
# DNDN_bar
# DNWT_bar

#convert the barplot to rose diagram
rose_ALL_bar<-ALL_bar + scale_y_continuous(breaks = c(0,10,20,30,40,50,60)) + coord_polar() + labs(x = "", y = "") + ggtitle("") +facet_wrap(~ treatment, nrow = 2) 
png("Rose plot, all.png")
plot(rose_ALL_bar)
dev.off()

#plot each treatment individually (Figure 2D)
#D18N_bar + scale_y_continuous(breaks=c(0,10,20),limits = c(0,20)) + coord_polar() + labs(x = "", y = "") + ggtitle("") + theme_bw() + theme(legend.position="none",panel.background = element_rect(fill="white",colour = "black", size = 0.75), axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.grid.major = element_line(colour = "gray",size=0.75), axis.title.x=element_blank())

#WT_bar + scale_y_continuous(breaks=c(0,10,20),limits = c(0,20)) + coord_polar() + labs(x = "", y = "") + ggtitle("") + theme_bw() + theme(legend.position="none",panel.background = element_rect(fill="white",colour = "black", size = 0.75), axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),panel.grid.major = element_line(colour = "gray",size=0.75),axis.title.x=element_blank())




##ANOVA to compare the effects of various infections on exPression of surface markers on CD4

# #Cluster 1
# AOV_C1<- aov (Cluster01~treatment, data = result_of_all)
# summary(AOV_C1)
# TukeyHSD(AOV_C1)
# 
# #Cluster 2
# AOV_C2<- aov (Cluster02~treatment, data = result_of_all)
# summary(AOV_C2)
# TukeyHSD(AOV_C2)
# 
# #Cluster 3
# AOV_C3<- aov (Cluster03~treatment, data = result_of_all)
# summary(AOV_C3)
# TukeyHSD(AOV_C3)
# 
# #Cluster 4
# AOV_C4<- aov (Cluster04~treatment, data = result_of_all)
# summary(AOV_C4)
# TukeyHSD(AOV_C4)
# 
# #Cluster 5
# AOV_C5<- aov (Cluster05~treatment, data = result_of_all)
# summary(AOV_C5)
# TukeyHSD(AOV_C5)
# 
# #Cluster 6
# AOV_C6<- aov (Cluster06~treatment, data = result_of_all)
# summary(AOV_C6)
# TukeyHSD(AOV_C6)
