#This file does cross-validation. Using data 09-14 to build up hclust, measuring block
setwd('/Users/lehathu/Desktop/AHC')
(list = ls())
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

#LOAD FUNCTION
source("fun.R")
## SIDE DATA
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]
short = read.csv("Public/NPI_zip.csv", colClasses = "character")

##IMPORT HRR DATA:
hrr = read.csv("Public/ZipHsaHrr14.csv", colClasses = 'character')
hrr = hrr[, c("zipcode14", "hrrnum")]
colnames(hrr) = c("zip", "hrr")
hrr = hrr[hrr$zip %in% zipcode$zip,]
hrr$hrr = as.numeric(paste(hrr$hrr))
hrr = merge(hrr, zipcode, by = 'zip', all.x =T)


load(file = "Data/ZipGraph_2014.RData")
B = calculateBlock(gzip, hrr, 'hrr')
sum(diag(B))/sum(B)


indeg = diag(B)/(colSums(B)+rowSums(B)-diag(B))
out = colnames(B)[which(indeg < 0.3)]
hrr[hrr$hrr %in% out,]



#DATA CREATED FROM "MakeNPIGraph.R"
load(file = "Data/HClust_0914.RData")
load(file = "Data/M_0914.RData") 


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



# GROUP NPI BY ZIP CODE
agg = aggregate(dt$clus, by = list(dt$zip), FUN = getmode)
colnames(agg) = c("zip", "clus")

agg = reorganize(agg, "clus") ## label of agg needs from 1:k
agg = merge(agg, zipcode, by = "zip", all.x = T)

length(unique(agg$clus))

load(file = 'Data/ZipGraph_0914.RData')
B = calculateBlock(gzip, agg, 'clus')
sum(diag(B))/sum(B)
indeg = diag(B)/(rowSums(B)+colSums(B)-diag(B))

## GREEDY
load(file = 'Data/ZipGraph_0914.RData')
g = induced_subgraph(gzip, V(gzip)[V(gzip)$name %in% agg$zip])
A = get.adjacency(g, attr = 'weight', sparse = T)
A = A + t(A)
A = A[order(colnames(A)), order(rownames(A))]
agg = agg[order(agg$zip),]
## greedy step
Z = makeMemMatrix(agg$clus)
Z <- as(Z, "dgCMatrix")
mem = Z
for(i in 1:1000){
	B = A %*% mem
	samp = sample(nrow(B), floor(nrow(B)/100), replace = F) # take 1% of doctors
	tmp = B[samp,]	
	col = apply(tmp, 1, function(z) reOrder(z))
	col = t(col)
	mem[samp,] = col
}
greedy = makeMemList1(mem)
agg = data.frame(agg,greedy) 
length(unique(greedy))


load('Data/ZipGraph_0914.RData')
B = calculateBlock(gzip, agg, 'greedy')
sum(diag(B)/sum(B))
indeg = diag(B)/colSums(B)

out = which(indeg < 0.4)
for(i in 1:length(out)){
	tmp = agg[agg$greedy == out[i],c('state', 'city')]
	print(unique(tmp))
}
#HRR
load('Data/ZipGraph_2014.RData')
B = calculateBlock(B, hrr, 'hrr')
sum(diag(B)/sum(B))


##USE KNN TO CLUSTER THE REST OF THE ZIPCODE:
knn = read.csv('Public/KnnZip_0914.csv', colClasses = 'character')
dat = agg[,c('zip','greedy', 'hrr' )]

knn = merge(knn, dat, all.x = T, by.x = 'zip1', by.y = 'zip' )
knn$greedy1 = knn$greedy
knn$greedy = NULL

dat = agg[,c('zip','greedy')]

knn = merge(knn, dat, all.x = T, by.x = 'zip2', by.y = 'zip' )
knn$greedy2 = knn$greedy
knn$greedy = NULL

knn = merge(knn, dat, all.x = T, by.x = 'zip3', by.y = 'zip' )
knn$greedy3 = knn$greedy
knn$greedy = NULL

knn = merge(knn, dat, all.x = T, by.x = 'zip4', by.y = 'zip' )
knn$greedy4 = knn$greedy
knn$greedy = NULL

knn = merge(knn, dat, all.x = T, by.x = 'zip5', by.y = 'zip' )
knn$greedy5 = knn$greedy
knn$greedy = NULL


x = as.matrix(cbind(knn$greedy1,knn$greedy2, knn$greedy3, knn$greedy4, knn$greedy5 ))

x = apply(x, 1, function(z) getmode(z))

knn = data.frame(knn$zip, knn$greedy1, knn$hrr)
colnames(knn) = c('zip', 'greedy', 'hrr')

dat = agg[,c('zip','greedy', 'hrr' )]
knn = rbind(knn, dat)
knn = merge(knn, zipcode, by = 'zip', all.x = T)

tab = as.data.frame(table(knn$greedy))
tab = tab[order(tab$Freq),]
head(tab)

#MOSAIC
knn1 = knn[!is.na(knn$hrr),]
knn1 = knn1[knn1$city == 'New York',]
knn1 = reorganize(knn1, 'greedy')
knn1 = reorganize(knn1, 'hrr')
loc = knn1[,c("longitude", "latitude")]
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = knn1$greedy
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

plot.dl(dl,clust)
outline <- map("usa", plot=FALSE) # returns a list of x/y coords
xbox = c(-1000,1000); ybox = c(-1000,1000)

polypath(c(outline$x, NA, c(xbox, rev(xbox))),
        c(outline$y, NA, rep(ybox, each=2)),
        col="white", rule="evenodd")
#map("county", add=T, lwd = .5)
legend('bottomright', legend=leglabs(round(brks,2)), fill=colors, bty="n",x.intersp = .5, y.intersp = .5, cex = 12, pt.cex = 2)

#See how NYC got chopped off
agg1 = agg[!is.na(agg$hrr),]
nyc = agg1[agg1$city == 'New York',]
table(nyc$hrr) #hrr groups all nyc together into 1 hrr
table(nyc$greedy)

agg1[agg1$greedy ==52, ] ## upper - weschester county
agg1[agg1$greedy ==79, ] ## bronx + part of long island
agg1[agg1$greedy ==94, ] ## manhattan downtown
agg1[agg1$greedy ==97, ] # queen
agg1[agg1$greedy ==98, ] #3 brooklyn
agg1[agg1$greedy ==103, ] 

#See how Madison got chopped off
mad = agg1[agg1$city == 'Madison' & agg1$state == 'WI',]
table(mad$hrr) #hrr groups all nyc together into 1 hrr
table(mad$greedy)

agg1[agg1$hrr == 299,]
agg1[agg1$greedy == 26,]
dt[dt$zip == '53711',]


spec = M[colnames(M) == 1740370519,]
spec = colnames(M)[spec>0]
spec = short[short$NPI %in% spec,]