##Combine all the doctor in the USA from 09-14 into graph. 

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
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

short = read.csv("Public/NPI_zip.csv", colClasses = "character")


makeGraphFromText = function(ref_filename, short){
	ref = read_csv(file = ref_filename, col_names = F,col_types = "cciii")
	colnames(ref) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")
	ref = ref[ref$NPI %in% short$NPI & ref$NPI2 %in% short$NPI,]
	g = graph.edgelist(as.matrix(ref[,1:2]), directed = F)
	E(g)$weight = ref$Unique
	return(g)
}

file_list = c("Public/physician-shared-patient-patterns-2010-days180.txt",
			"Public/physician-shared-patient-patterns-2011-days180.txt",
			"Public/physician-shared-patient-patterns-2012-days180.txt",
			"Public/physician-shared-patient-patterns-2013-days180.txt",
			"Public/physician-shared-patient-patterns-2014-days180.txt")

g = makeGraphFromText("Public/physician-shared-patient-patterns-2009-days180.txt", short)

for(i in 1:length(file_list)){
	gtemp = makeGraphFromText(file_list[i], short)
	g = g %u% gtemp
	
	E(g)$weight_1[is.na(E(g)$weight_1)] = 0
	E(g)$weight_2[is.na(E(g)$weight_2)] = 0
	
	E(g)$weight = E(g)$weight_1 + E(g)$weight_2
	
	g <- remove.edge.attribute(g, "weight_1")
	g <- remove.edge.attribute(g, "weight_2")
	
	print(length(V(gtemp)))
	print(i)
}

W = get.adjacency(g, attr = 'weight', sparse = T)

# Regularize W
O = rowSums(W)
tau = rep(mean(O), length(O))

L_reg = Diagonal(length(O), 1/sqrt(O+tau))%*%W%*%Diagonal(length(O), 1/sqrt(O+tau))

gref = graph.adjacency(L_reg, weighted = T, diag = F)
core = coreness(gref)
hist(core)
gsm = induced_subgraph(gref,core>median(core))
M = get.adjacency(gsm,attr = "weight", sparse = T)
save(M, file = "Data/M_0914.RData")

##SpareAHC
hclust = sparseAHC(M, linkage = "average")
save(hclust, file = "Data/HClust_0914.RData")

