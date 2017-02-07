#This file clusters usa into hsa (over 3k)
setwd('/Users/lehathu/Desktop/AHC')
rm(list = ls())
library(Matrix)
library(sparseAHC)
library(zipcode)
data(zipcode)
library(geosphere)
library(igraph)
library(deldir)
library(RColorBrewer)
library(maps)
library(classInt)
library(maptools)
library(rARPACK)
library(Matrix)
##CUSTOM FUNCTION

# reorganize the clus
reorganize = function(df, clus){
	col = which(colnames(df) == clus)
	unique_label = sort(unique(df[,col]))
	ind = c(1:length(unique(df[,col])))
	mapping = data.frame(unique_label,ind)

	for(i in 1:length(mapping$unique_label)){
		df[,col][df[,col] == mapping$unique_label[i]] = mapping$ind[i]
	}
	return(df)
}

#function to group NPI by zipcode
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

#greedy-step functions
makeMemMatrix = function(clus){
	x = sort(unique(clus))
	m = matrix(0, nrow = length(clus), ncol = length(x))
	for (i in 1:length(x)){
		m[,i][clus==x[i]] = 1
	}
	return(m)
}
makeMemList = function(m){
	clus = apply(m, 1, function(x) which(x>0))
	return(clus)
}

reOrder = function(vect){
	a = rep(0, length(vect))
	a[which.max(vect)] = 1
	return(a)
}

calculateBlockHsa = function(g, df, clus){
	col = which(colnames(df) == clus)
	W = get.adjacency(g, attr = 'weight', sparse = T)
	W = W[rownames(W) %in% df$hsa, colnames(W) %in% df$hsa]
	W = W[order(rownames(W)),order(colnames(W))]
	
	df = df[df$hsa %in% colnames(W),]
	df = df[order(df$hsa),]
	
	Z = makeMemMatrix(df[,col])
	B = t(Z) %*% W %*% Z
	B = as.matrix(B)
	return(B)
}

calculateBlock = function(graph, df, clus){
	col = which(colnames(df) == clus)
	W = get.adjacency(graph, attr = 'weight', sparse = T)
	W = W[rownames(W) %in% df$zip, colnames(W) %in% df$zip]
	W = W[order(rownames(W)),order(colnames(W))]
	
	df = df[df$zip %in% colnames(W),]
	df = df[order(df$zip),]
	
	Z = makeMemMatrix(df[,col])
	B = t(Z) %*% W %*% Z
	B = as.matrix(B)
	colnames(B) = sort(unique(df[,col]))
	return(B)
}




## SIDE DATA
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]
short = read.csv("Public/NPI_zip.csv", colClasses = "character")

##IMPORT HRR DATA:
hrr = read.csv("Public/ZipHsaHrr14.csv", colClasses = 'character')
colnames(hrr) = c("zip", "hsa", 'hsacity', 'hsastate', 'hrr', 'hrrcity', 'hrrstate')s
hrr = hrr[hrr$zip %in% zipcode$zip,]

hrr$hrr = as.numeric(paste(hrr$hrr))
hrr$hsa = as.numeric(paste(hrr$hsa))
hrr = merge(hrr, zipcode, by = 'zip', all.x =T)
length(unique(hrr$hsa))

#shrink hrr
hsa = hrr[, c('hsa', 'hrr')]
hsa = unique(hsa)

load(file = "Data/hsaGraph.RData")
load(file = "Data/M_hsa.RData")
load(file = "Data/HClust_hsa.RData")

## see how well hrr cut hsa
B = calculateBlock(g, hrr, 'hrr')

##CUT THE DENDROGRAM AND ASSIGN NPI 

k = 340
x = cutree(hclust,k)

dt = data.frame(colnames(M), x)
colnames(dt) = c("NPI", "clus")
dt = merge(dt, short, by = "NPI", all.x = T)
dt = merge(dt, zipcode, by = "zip", all.x = T)

tab = as.data.frame(table(dt$clus))
tab = tab[order(tab$Freq),]
head(tab)
