rm(list = ls())
library(Matrix)
library(sparseAHC)
library(zipcode)
data(zipcode)
library(geosphere)

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
short = read.csv("Data/NPI_zip.csv", colClasses = "character")

#DATA CREATED FROM "MakeNPIGraph.R"
load(file = "Data/HClust_NPI.RData")
load(file = "Data/NPI_Matrix.RData")
NPI = read.csv("Data/HClust_NPI.csv")


##cluster into HRR
k = 306
x = cutree(hclust,k)

## greedy step
Z = makeMemMatrix(x)
Z <- as(Z, "dgCMatrix")
mem = Z
for(i in 1:10000){
	B = M %*% mem
	samp = sample(nrow(B), floor(nrow(B)/1000), replace = F) # take 1% of doctors
	tmp = B[samp,]	
	col = apply(tmp, 1, function(z) reOrder(z))
	col = t(col)
	mem[samp,] = col
}
clus = makeMemList(mem)
sum(clus!=x)

# group the NPI by Zipcode
dt = data.frame(NPI, x)
colnames(dt) = c("NPI", "clus")
dt = merge(dt, short, by = "NPI", all.x = T)
dt = merge(dt, zipcode, by = "zip", all.x = T)
agg = aggregate(dt$clus, by = list(dt$zip), FUN = getmode)
colnames(agg) = c("zip", "clus")
length(unique(agg$clus))

agg = reorganize(agg, "clus") ## label of agg needs from 1:k
agg = merge(agg, zipcode, by = "zip", all.x = T)

print("Process started")
#Assigning zipcode to cluster
excl = zipcode[!(zipcode$zip %in% agg$zip), ]
excl$clus = 0
for(i in 1:dim(excl)[1]){
	tmp = excl[i,]
	loc = c(tmp$longitude, tmp$latitude)
    m = apply(cbind(agg$longitude, agg$latitude), 1, function(p, loc) distCosine(p,loc), loc = loc)
	excl$clus[i] = agg$clus[which.min(m)]
	if(i%%10==0){print(i)}
}
excl = excl[, colnames(agg)]
write.csv(excl, "Data/AllZip_clus.csv", row.names = F)
write.csv(agg, "Data/CoreZip_clus.csv", row.names = F)

agg = read.csv("Data/CoreZip_clus.csv")
excl = read.csv("Data/AllZip_clus.csv")
agg = rbind(agg, excl)

##Mosaic plot
loc = agg[,c("longitude", "latitude")]
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = agg$clus
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/AllZip_HRR", max(agg$clus), "_greedy.pdf", sep = "")
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

hrr = read.csv("Data/ZipHsaHrr14.csv")
hrr = hrr[, c("zipcode14", "hrrnum")]
colnames(hrr) = c("zip", "hrr")

agg = merge(agg, hrr, by = "zip", all.x = T)
agg = agg[!is.na(agg$hrr),]
compare(agg$clus, agg$hrr, method = "rand")
compare(agg$clus, agg$hrr, method = "nmi")
compare(agg$clus, agg$hrr, method = "adjusted.rand")
