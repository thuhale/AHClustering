##Build up zip graph from referral data in 09-14
rm(list = ls())
library(readr)
library(igraph)
library(Matrix)
library(sparseAHC)
library(MCL)
library(zipcode)
data(zipcode)
library(maps)
#setwd('/Users/lehathu/Desktop/AHC')

short = read.csv("Public/NPI_zip.csv", colClasses = "character")


makeZipGraph = function(ref_filename, short){
	ref = read_csv(file = ref_filename, col_names = F,col_types = "cciii")
	colnames(ref) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")
	ref = ref[ref$NPI %in% short$NPI & ref$NPI2 %in% short$NPI,]
	
	ref = merge(ref, short, by = "NPI", all.x = T)
	ref = merge(ref, short, by.x = "NPI2", by.y = "NPI", all.x = T)
	colnames(ref)[6:7] = c("zip", "zip2")
	ref = ref[, c(2,1,3,4,5,6,7)]
	ref$zip = clean.zipcodes(ref$zip)
	ref$zip2 = clean.zipcodes(ref$zip2)

	#make the graph
	gzip = graph.edgelist(as.matrix(ref[,6:7]), directed = F)
	E(gzip)$weight = ref$Unique
	gzip = simplify(gzip, remove.loops = FALSE)
	
	return(gzip)
}

file_list = c("Public/physician-shared-patient-patterns-2010-days180.txt",
			"Public/physician-shared-patient-patterns-2011-days180.txt",
			"Public/physician-shared-patient-patterns-2012-days180.txt",
			"Public/physician-shared-patient-patterns-2013-days180.txt",
			"Public/physician-shared-patient-patterns-2014-days180.txt")

gzip = makeZipGraph("Public/physician-shared-patient-patterns-2009-days180.txt", short)

for(i in 1:length(file_list)){
	gtemp = makeZipGraph(file_list[i], short)
	gzip = gzip %u% gtemp
	
	E(gzip)$weight_1[is.na(E(gzip)$weight_1)] = 0
	E(gzip)$weight_2[is.na(E(gzip)$weight_2)] = 0
	
	E(gzip)$weight = E(gzip)$weight_1 + E(gzip)$weight_2
	
	gzip <- remove.edge.attribute(gzip, "weight_1")
	gzip <- remove.edge.attribute(gzip, "weight_2")
	
}
save(gzip, file = "Data/ZipGraph_0914.RData")


