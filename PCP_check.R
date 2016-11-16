rm(list = ls())
library(igraph)
setwd("/Users/lehathu/Desktop/AHC")
library(Matrix)
library(sparseAHC)
library(MCL)
library(zipcode)
data(zipcode)
library(maps)
library(deldir)
library(RColorBrewer)
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-US continential 
zipcode = zipcode[!(zipcode$state %in% peri),]
## This file reads in PCP edge_list network "PCP_Referral_List.csv" generated from MakePCPGraph.R file and performing cluster on that. 
## It also plots the cluster on map. Doctor's location is found on "NPI_zip.csv"
# read in doctor's zipcode
short = read.csv("Data/NPI_zip.csv", colClasses = "character")
#create a merge matrix
load("Data/PCP_Referral_graph.RData")
M = get.adjacency(gsm,attr = "weight", sparse = T)

x = sample(dim(M)[1], dim(M)[1], replace = F)
Mb = M[x,x]

hclust = sparseAHC(M, linkage = "average")
hclust_b = sparseAHC(Mb, linkage = "average")
##Shuffling to see if AHC is stable
k = 738
k = 256
x = cutree(hclust,k)
x_b = cutree(hclust_b,k)

tab1 = as.data.frame(table(x))
tab2 = as.data.frame(table(x_b))

small1 = tab1$x[tab1$Freq == 1]
small2 = tab2$x[tab2$Freq == 1]

nei = neighbors(gsm,small1[4])
table(x[colnames(M) %in% nei$name])

sort(colnames(M)[x %in% small1])
sort(colnames(Mb)[x_b %in% small2])

##check the degree of those individual NPI
sort(which(x %in% small1))
v = 36051
degree(gsm, v)
ego = induced_subgraph(gsm, neighborhood(gsm,2,373)[[1]])
heatmap(as.matrix(get.adjacency(ego)))


hist(degree(gsm, which(x %in% small1))

which.min(degree(gsm))
samp = sample(length(x), 10000)
xx = x[samp]
yy = degree(gsm)[samp]

#degree
plot(xx, yy)
lines(lowess(xx,yy), col = "red", lwd =4)
plot(xx, log(yy))
lines(lowess(xx,log(yy)), col = "red", lwd =4)

#coreness
yy =coreness(gsm)[samp]
plot(xx, log(yy))
lines(lowess(xx,log(yy)), col = "red", lwd =4)

#Transitivity
tmp = transitivity(gsm, type = "local")
yy = tmp[samp]
plot(xx, yy)
lines(lowess(xx,yy), col = "red", lwd =4)

#Betwennness
tmp = betweenness(gsm, V(gsm)[samp], directed = F)