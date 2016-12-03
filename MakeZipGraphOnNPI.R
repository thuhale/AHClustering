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


agg = read.csv("Data/CoreZip_clus.csv", colClasses = "character")
hrr = read.csv("Data/ZipHsaHrr14.csv")
hrr = hrr[, c("zipcode14", "hrrnum")]
colnames(hrr) = c("zip", "hrr")
hrr$zip = clean.zipcodes(hrr$zip)

agg = merge(agg, hrr, by = "zip", all.x = T)
#agg = agg[!is.na(agg$hrr),]
agg = reorganize(agg, "hrr") ## label of agg needs from 1:k

##MAKE ZIP GRAPH
short = read.csv("Data/NPI_zip.csv")

ref = read_csv(file = "Data/physician-shared-patient-patterns-2014-days180.txt", col_names = F,col_types = "cciii")

colnames(ref) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")
ref = ref[ref$NPI %in% short$NPI & ref$NPI2 %in% short$NPI,]

ref = merge(ref, short, by = "NPI", all.x = T)
ref = merge(ref, short, by.x = "NPI2", by.y = "NPI", all.x = T)
colnames(ref)[6:7] = c("zip", "zip2")
ref = ref[, c(2,1,3,4,5,6,7)]
ref$zip = clean.zipcodes(ref$zip)
ref$zip2 = clean.zipcodes(ref$zip2)

tmp = ref[ref$zip %in% agg$zip & ref$zip2 %in% agg$zip, ]
ref = tmp
#make the graph
g = graph.edgelist(as.matrix(ref[,6:7]), directed = F)
E(g)$weight = ref$Unique
g = simplify(g)
save(g, file = "Data/Zip_Cluster.RData")