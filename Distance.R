#This file does cross-validation. Using data 09-14 to build up hclust, measuring block
setwd('/Users/lehathu/Desktop/AHC')
library(vegan)
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
source("functions.R")
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





#DATA CREATED FROM "MakeNPIGraph.R"
load(file = "Data/HClust_0914.RData")
load(file = "Data/M_0914.RData") 


##CUT THE DENDROGRAM AND ASSIGN NPI 

k = 338
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

load(file = 'Data/ZipGraph_2014.RData')
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
agg = reorganize(agg, 'greedy')

agg$col = as.numeric(substr(agg$zip, 1,1))
col = agg[,c('greedy', 'col')]
col = aggregate(agg$col, by = list(agg$greedy), FUN = getmode)
colnames(col) = c('greedy', 'col')

load(file = 'Data/ZipGraph_2014.RData')
B = calculateBlock(gzip, agg, 'greedy')
sum(diag(B))/sum(B)


pdat = mdsPlot(agg, 'greedy', 'zip', gzip, 'longitude', 'latitude')
pdat_hrr = mdsPlot(hrr, 'hrr', 'zip', gzip, 'longitude', 'latitude')



coord1 = as.matrix(pdat[, c('x', 'y')])
coord1[,2] = -1*coord1[,2]

coord2 = as.matrix(pdat[, c('long', 'lat')])


#coord2[,2] = -1*coord2[,2]

sub = procrustes(coord1, coord2, scale = TRUE, symmetric = FALSE, scores = "sites")
ob = plot(sub, kind = 1)

