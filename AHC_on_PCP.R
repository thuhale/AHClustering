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
#Point-wise plot
#k = 738
k = 256
x = cutree(hclust,k)
x_b = cutree(hclust_b,k)


pdf(file = paste("Pix/PCP_Hist", max(x), ".pdf", sep = ""))
hist(x, main = paste("PCP_Hist", max(x)))
dev.off()

dt = data.frame(colnames(M), x)
colnames(dt) = c("NPI", "clus")
dt = merge(dt, short, by = "NPI", all.x = T)
dt = merge(dt, zipcode, by = "zip", all.x = T)

#group dt by zipcode
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

agg = aggregate(dt$clus, by = list(dt$zip), FUN = getmode)
colnames(agg) = c("zip", "clus")
length(unique(agg$clus))

# reorganize the clus
unique_label = sort(unique(agg$clus))
ind = c(1:length(unique(agg$clus)))
mapping = data.frame(unique_label,ind)

for(i in 1:length(mapping$unique_label)){
	agg$clus[agg$clus == mapping$unique_label[i]] = mapping$ind[i]
}

agg = merge(agg, zipcode, by = "zip", all.x = T)
loc =agg[,c("longitude", "latitude")]

##Mosaic plot
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = agg$clus
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/PCP_HRR", max(dt$clus), ".pdf", sep = "")
pdf(file = filename, height = 50, width = 100)
plot.dl(dl,clust)
outline <- map("usa", plot=FALSE) # returns a list of x/y coords
xbox = c(-1000,1000); ybox = c(-1000,1000)
 # par("usr")[1:2]; ybox = par("usr")[3:4]
# create the grid path in the current device
polypath(c(outline$x, NA, c(xbox, rev(xbox))),
        c(outline$y, NA, rep(ybox, each=2)),
        col="white", rule="evenodd")
map("county", add=T, lwd = .5)
dev.off()

##Data export to Python
order = order.dendrogram(as.dendrogram(hclust))
order = sort(order, index.return = T)$ix
sorted_order = data.frame(colnames(M), 1:length(colnames(M)), order)
colnames(sorted_order) = c("NPI", "index", "order")
sorted_order = merge(sorted_order, short, by = "NPI", all.x = T)
sorted_order = merge(sorted_order, zipcode, by = "zip", all.x = T)
sorted_order = sorted_order[,1:6]
sorted_order = sorted_order[order(sorted_order$index),]
merge = hclust$merge
merge = as.data.frame(merge)
colnames(merge) = c("V1", "V2")

write.csv(sorted_order, "Output/pcp_sorted_order.csv", row.names = F)
write.csv(merge, "Output/pcp_merge.csv", row.names = F)
#Feed back from Python
x = read.csv("pcp_label.csv", header = F)
dt = data.frame(colnames(M), x$V2)
colnames(dt) = c("NPI", "clus")
dt = merge(dt, short, by = "NPI", all.x = T)
dt = merge(dt, zipcode, by = "zip", all.x = T)

pdf(file = paste("Pix/PCP_Hist", max(dt$clus), "_balanced.pdf",".pdf",sep = ""))
hist(dt$clus, main = paste("PCP_Hist", max(dt$clus), "Balance"))
dev.off()
#group dt by zipcode
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

agg = aggregate(dt$clus, by = list(dt$zip), FUN = getmode)
colnames(agg) = c("zip", "clus")
length(unique(agg$clus))

# reorganize the clus
unique_label = sort(unique(agg$clus))
ind = c(1:length(unique(agg$clus)))
mapping = data.frame(unique_label,ind)

for(i in 1:length(mapping$unique_label)){
	agg$clus[agg$clus == mapping$unique_label[i]] = mapping$ind[i]
}

agg = merge(agg, zipcode, by = "zip", all.x = T)
##Mosaic plot
loc = dt[,c("longitude", "latitude")]
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = dt$clus
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/PCP_HRR", max(dt$clus), "_balanced.pdf", sep = "")
pdf(file = filename, height = 50, width = 100)
plot.dl(dl,clust)
outline <- map("usa", plot=FALSE) # returns a list of x/y coords
xbox = c(-1000,1000); ybox = c(-1000,1000)
 # par("usr")[1:2]; ybox = par("usr")[3:4]
# create the grid path in the current device
polypath(c(outline$x, NA, c(xbox, rev(xbox))),
        c(outline$y, NA, rep(ybox, each=2)),
        col="white", rule="evenodd")
map("county", add=T, lwd = .5)
dev.off()

map("state", col = 1)
for(i in 1:length(unique(dt$clus))){
	points(dt$longitude[dt$clus==i], dt$latitude[dt$clus==i], col = i, pch = 19, cex = 0.5)
}