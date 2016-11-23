rm(list = ls())
library(Matrix)
library(sparseAHC)
library(MCL)
library(zipcode)
data(zipcode)
library(maps)
library(deldir)
library(igraph)
library(RColorBrewer)


peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

short = read.csv("Data/NPI_zip.csv", colClasses = "character")

load("Data/NPI_Matrix.RData")
load("Data/HClust_NPI.RData") #load hclust
NPI = read.csv("Data/HClust_NPI.csv")

k = 306
x = cutree(hclust,k)

Z = makeMemMatrix(x)
Z <- as(Z, "dgCMatrix")
mem = Z
for(i in 1:1000){
	B = M %*% mem
	samp = sample(nrow(B), floor(nrow(B)/100), replace = F)
	tmp = B[samp,]	
	col = apply(tmp, 1, function(z) reOrder(z))
	col = t(col)
	mem[samp,] = col
}
clus = makeMemList(mem)
# group the NPI by Zipcode
dt = data.frame(NPI, x)
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
agg = reorganize(agg, "clus")
agg = merge(agg, zipcode, by = "zip", all.x = T)

##Mosaic plot
loc = agg[,c("longitude", "latitude")]
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = agg$clus
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/PCP_HRR", max(dt$clus), "_greedy.pdf", sep = "")
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
