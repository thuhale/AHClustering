rm(list = ls())
library(readr)
library(Matrix)
library(sparseAHC)
library(zipcode)
data(zipcode)
library(geosphere)
library(igraph)
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

##cluster into HRR
k = 306
x = cutree(hclust,k)

## greedy step
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


tab = as.data.frame(table(dt$clus))
tab = tab[order(tab$Freq),]

tab2 = as.data.frame(table(dt$greedy))
tab2 = tab2[order(tab2$Freq),]


# group the NPI by Zipcode
agg = aggregate(dt$clus, by = list(dt$zip), FUN = getmode)
colnames(agg) = c("zip", "clus")

agg = reorganize(agg, "clus") ## label of agg needs from 1:k
agg = merge(agg, zipcode, by = "zip", all.x = T)

length(unique(agg$clus))
tab = as.data.frame(table(agg$clus))
tab = tab[order(tab$Freq),]
head(tab)
W = get.adjacency(gzip, attr = 'weight', sparse = T)
W = W[rownames(W) %in% agg$zip, colnames(W) %in% agg$zip]
W = W[order(rownames(W)),order(colnames(W))]
agg = agg[order(agg$zip),]
A = W
diag(A) = 0
## greedy step
Z = makeMemMatrix(agg$clus)
Z <- as(Z, "dgCMatrix")
mem = Z
nclu = c()
for(i in 1:10000){
	B = A %*% mem
	samp = sample(nrow(B), floor(nrow(B)/100), replace = F) # take 1% of doctors
	tmp = B[samp,]	
	col = apply(tmp, 1, function(z) reOrder(z))
	col = t(col)
	mem[samp,] = col
	nclu = c(nclu,sum(colSums(mem)!=0))
}
greedy = makeMemList(mem)
agg = data.frame(agg, greedy)
length(unique(agg$greedy))

tab3 = as.data.frame(table(agg$greedy))
tab3 = tab3[order(tab3$Freq),]
head(tab3)

Z = makeMemMatrix(agg$greedy)
B = t(Z) %*% W %*% Z 
B = as.matrix(B)
heatmap(B)
#Number of within cluster compared to outside cluster
sum(diag(B))/sum(B)
diag(B)/rowSums(B)
hist(diag(B)/rowSums(B))

agg[agg$clus == which.min((diag(B)/rowSums(B))),]
agg[agg$clus == which.max((diag(B)/rowSums(B))),]

## Compare it with hrr cluster
hrr = read.csv("Public/ZipHsaHrr14.csv")
hrr = hrr[, c("zipcode14", "hrrnum")]
colnames(hrr) = c("zip", "hrr")
hrr$zip = clean.zipcodes(hrr$zip)

agg = merge(agg, hrr, by = "zip", all.x = T)
agg1 = agg[!is.na(agg$hrr),]

agg1 = reorganize(agg1, "hrr")
Z1 = makeMemMatrix(agg1$hrr)
W1 = W[rownames(W) %in% agg1$zip, colnames(W) %in% agg1$zip]
W1 = W1[order(rownames(W1)),order(colnames(W1))]
agg1 = agg1[order(agg1$zip),]

length(unique(agg1$hrr))
tab1 = as.data.frame(table(agg1$hrr))
tab1 = tab1[order(tab1$Freq),]
head(tab1)

B1 = t(Z1) %*% W1 %*% Z1 
B1 = as.matrix(B1)
heatmap(B1)
#Number of within cluster compared to outside cluster
sum(diag(B1))/sum(B1)
summary(diag(B1)/rowSums(B1))
hist(diag(B)/rowSums(B))
agg1[agg1$hrr == which.min((diag(B1)/rowSums(B1))),]



#Random Walk
Z = makeMemMatrix(agg$clus)
B = t(Z) %*% W %*% Z
O = rowSums(B)
B = diag(1/O) %*% B
B = as.matrix(B)
pdf(file = "Pix/BlockMatrix_RW.pdf")
heatmap(B)
dev.off()

self = diag(B)
tmp = B
diag(tmp) = 0
O = rowSums(tmp)

x = self - O
agg[agg$clus == which.max(x),]
agg[agg$hrr==,]

agg[agg$city == "Detroit" & agg$state == "MI",]