rm(list = ls())
library(data.table)
library(readr)
library(igraph)
library(Matrix)
library(sparseAHC)
library(MCL)
library(zipcode)
data(zipcode)
library(maps)
library(maptools)
library(geosphere)
library(mclust)
library(mcclust)

peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

short = read.csv("Data/NPI_zip.csv", colClasses = "character")

load("Data/NPI_Matrix.RData")
load("Data/HClust_NPI.RData") #load hclust
NPI = read.csv("Data/HClust_NPI.csv")

k = 306
x = cutree(hclust,k)
hist(x)

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

excl = read.csv("Data/AllZip_clus.csv", colClasses = "character")
dim(excl)
excl = excl[, colnames(agg)]
excl$clus = as.numeric(paste(excl$clus))
excl$latitude = as.numeric(paste(excl$latitude))
excl$longitude = as.numeric(paste(excl$longitude))
excl = excl[, colnames(agg)]
excl = rbind(agg, excl)

#Mosaic plot
loc =excl[,c("longitude", "latitude")]
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = excl$clus
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/NPI_HRR_USA_AllZipp", max(excl$clus), ".pdf", sep = "")
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

hrr= read.csv("Data/ZipHsaHrr14.csv", colClasses ="character")
hrr = hrr[, c(1,5)]
colnames(hrr) = c("zip", "hrr")
hrr$hrr = as.numeric(paste(hrr$hrr))
hrr1 = hrr[hrr$zip %in% agg$zip,]
hrr1 = reorganize(hrr1, "hrr")

agg1 = agg[agg$zip %in% hrr1$zip,]
agg1 = merge(agg1, hrr1, by = "zip", all.x = T)
compare(agg1$clus, agg1$hrr, method = c("rand"))
compare(agg1$clus, agg1$hrr, method = c("adjusted.rand"))
compare(agg1$clus, agg1$hrr, method = c("nmi"))

hrr1 = hrr[hrr$zip %in% excl$zip,]
hrr1 = reorganize(hrr1, "hrr")

excl1 = excl[excl$zip %in% hrr1$zip,]
excl1 = merge(excl1, hrr1, by = "zip", all.x = T)
compare(excl1$clus, excl1$hrr, method = c("rand"))
compare(excl1$clus, excl1$hrr, method = c("adjusted.rand"))
compare(excl1$clus, excl1$hrr, method = c("nmi"))


#Mosaic plot
loc =agg[,c("longitude", "latitude")]
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = agg$clus
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/HRR_by_Medicareshfiosdhfo", max(excl1$clus), ".pdf", sep = "")
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

