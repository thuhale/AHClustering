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
	m = matrix(0, nrow = length(clus), ncol = length(unique(clus)))
	for (i in 1:length(unique(clus))){
		m[,i][clus==i] = 1
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

## SIDE DATA
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]
short = read.csv("Public/NPI_zip.csv", colClasses = "character")

#DATA CREATED FROM "MakeNPIGraph.R"
load(file = "Data/HClust_2014.RData")
load(file = "Data/M_2014.RData") 
load(file = "Data/ZipGraph_2014.RData")

##CUT THE DENDROGRAM AND ASSIGN NPI 

k = 306
x = cutree(hclust,k)

# greedy step
Z = makeMemMatrix(x)
Z <- as(Z, "dgCMatrix")
mem = Z
for(i in 1:1000){
	B = M %*% mem
	samp = sample(nrow(B), floor(nrow(B)/100), replace = F) # take 1% of doctors
	tmp = B[samp,]	
	col = apply(tmp, 1, function(z) reOrder(z))
	col = t(col)
	mem[samp,] = col
}
clus = makeMemList(mem)
dt = data.frame(colnames(M), x, clus)
colnames(dt) = c("NPI", "clus", "greedy")
dt = merge(dt, short, by = "NPI", all.x = T)
dt = merge(dt, zipcode, by = "zip", all.x = T)


tab1 = as.data.frame(table(dt$clus))
tab1 = tab1[order(tab1$Freq),]
head(tab1)

tab2 = as.data.frame(table(dt$greedy))
tab2 = tab2[order(tab2$Freq),]
head(tab2)

# GROUP NPI BY ZIP CODE
agg = aggregate(dt$clus, by = list(dt$zip), FUN = getmode)
colnames(agg) = c("zip", "clus")

agg = reorganize(agg, "clus") ## label of agg needs from 1:k
agg = merge(agg, zipcode, by = "zip", all.x = T)

length(unique(agg$clus))
tab3 = as.data.frame(table(agg$clus))
tab3 = tab3[order(tab3$Freq),]
head(tab3)


#USE KNN TO CLUSTER THE REST OF THE ZIPCODE:
knn = read.csv('/Users/lehathu/Desktop/AHC/Public/KnnZip.csv')




