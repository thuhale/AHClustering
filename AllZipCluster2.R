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

k = 309
x = cutree(hclust,k)

dt = data.frame(colnames(M), x)
colnames(dt) = c("NPI", "clus")
dt = merge(dt, short, by = "NPI", all.x = T)
dt = merge(dt, zipcode, by = "zip", all.x = T)


tab1 = as.data.frame(table(dt$clus))
tab1 = tab1[order(tab1$Freq),]
head(tab1)


# GROUP NPI BY ZIP CODE
agg = aggregate(dt$clus, by = list(dt$zip), FUN = getmode)
colnames(agg) = c("zip", "clus")

agg = reorganize(agg, "clus") ## label of agg needs from 1:k
agg = merge(agg, zipcode, by = "zip", all.x = T)

length(unique(agg$clus))
tab3 = as.data.frame(table(agg$clus))
tab3 = tab3[order(tab3$Freq),]
head(tab3)


## GREEDY
g = induced_subgraph(gzip, V(gzip)[V(gzip)$name %in% agg$zip])
A = get.adjacency(g, attr = 'weight', sparse = T)
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
greedy = makeMemList(mem)
agg = data.frame(agg,greedy) 
length(unique(greedy))
agg = reorganize(agg, 'greedy')
tm = data.frame(agg$zip, agg$greedy)
colnames(tm) = c('zip', 'zip_clus')
dt = merge(dt, tm, by = 'zip', all.x = T)

tab4 = as.data.frame(table(dt$zip_clus))
tab4 = tab4[order(tab4$Freq),]
head(tab4)
hist(dt$zip_clus, main = 'Histogram of cluster size', ylab = 'Number of NPI', xlab = 'cluster')

##BLOCK MODEL:
W = get.adjacency(gzip, attr = 'weight', sparse = T)
W = W[rownames(W) %in% agg$zip, colnames(W) %in% agg$zip]
W = W[order(rownames(W)),order(colnames(W))]
agg = agg[order(agg$zip),]

Z = makeMemMatrix(agg$greedy)

B = t(Z) %*% W %*% Z 
B = as.matrix(B)
#png(file = 'Pix/block.png')
#heatmap(B)
#dev.off()
#sum(diag(B))/sum(rowSums(B))

#Number of within cluster compared to outside cluster
indeg = diag(B)/rowSums(B)
#summary(indeg)
#hist(indeg, main = "Histogram of the clusters's endogeneity", xlab = '% of in-clustered co-ocurrences')
#plot(indeg ~ colSums(Z))

# Explore clusters that are highly exogeneous
small = which(indeg < 0.4)
agg[agg$clus == small,]
dt[dt$zip_clus %in% small,]



#SIMILARITY MEASURE W/HRR
hrr = read.csv("Public/ZipHsaHrr14.csv", colClasses = 'character')
hrr = hrr[, c("zipcode14", "hrrnum")]
colnames(hrr) = c("zip", "hrr")

agg = merge(agg, hrr, by = "zip", all.x = T)
agg1 = agg[!is.na(agg$hrr),]

compare(agg1$greedy, agg1$hrr, method = "rand")
compare(agg1$greedy, agg1$hrr, method = "nmi")
compare(agg1$greedy, agg1$hrr, method = "adjusted.rand")

## hrr sizes 
dt = merge(dt, hrr, by = 'zip', all.z = T)
dt1 = dt[!is.na(dt$hrr),]
dt1$hrr = as.numeric(dt1$hrr)
dt1 = reorganize(dt1, 'hrr')

hist(dt1$hrr, main = 'Histogram of HRR size', ylab = 'Number of NPI', xlab = 'HRR')

#block model for hrr
W = W[rownames(W) %in% agg1$zip, colnames(W) %in% agg1$zip]
W = W[order(rownames(W)),order(colnames(W))]

agg1 = agg1[order(agg1$zip),]
agg1$hrr = as.numeric(agg1$hrr)
agg1 = reorganize(agg1, 'hrr')
Z = makeMemMatrix(agg1$hrr)

B = t(Z) %*% W %*% Z 
B = as.matrix(B)

sum(diag(B))/sum(B)
hrr_indeg = diag(B)/rowSums(B)

#hist(diag(B)/rowSums(B), main = "Histogram of the HRR's endogeneity", xlab = '% of in-HRR co-ocurrences')

##USE KNN TO CLUSTER THE REST OF THE ZIPCODE:
knn = read.csv('Public/KnnZip_2014.csv', colClasses = 'character')
dat = agg1[,c('zip','greedy', 'hrr' )]

knn = merge(knn, dat, all.x = T, by.x = 'zip1', by.y = 'zip' )
knn$greedy1 = knn$greedy
knn$greedy = NULL

dat = agg1[,c('zip','greedy')]

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

knn = data.frame(knn$zip, x, knn$hrr)
colnames(knn) = c('zip', 'greedy', 'hrr')

dat = agg1[,c('zip','greedy', 'hrr' )]

knn = rbind(knn, dat)

knn = merge(knn, zipcode, by = 'zip', all.x = T)
knn = knn[!is.na(knn$hrr),]


#MOSAIC PLOT
##Mosaic plot

brks<-classIntervals(indeg, n=9, style="quantile")
brks<- brks$brks
colors <- brewer.pal(9, "YlOrRd")

loc = knn[,c("longitude", "latitude")]
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = knn$hrr
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/AllZip_HRR_Constrast_Scaled", max(agg$greedy), '.pdf', sep = "")
pdf(file = filename, height = 50, width = 100)

plot.dl(dl,clust, indeg = hrr_indeg, ruler = indeg)
outline <- map("usa", plot=FALSE) # returns a list of x/y coords
xbox = c(-1000,1000); ybox = c(-1000,1000)
 # par("usr")[1:2]; ybox = par("usr")[3:4]
# create the grid path in the current device
polypath(c(outline$x, NA, c(xbox, rev(xbox))),
        c(outline$y, NA, rep(ybox, each=2)),
        col="white", rule="evenodd")
#map("county", add=T, lwd = .5)
legend('bottomright', legend=leglabs(round(brks,2)), fill=colors, bty="n",x.intersp = .5, y.intersp = .5, cex = 12, pt.cex = 2)
dev.off()


