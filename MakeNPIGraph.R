rm(list = ls())
library(data.table)
library(readr)
library(igraph)
setwd("/Users/lehathu/Desktop/AHC")
library(Matrix)
library(sparseAHC)
library(MCL)
library(dplyr)
library(zipcode)
data(zipcode)
library(maps)
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

short = read.csv("Data/NPI_zip.csv")
ref = read_csv(file = "Data/physician-shared-patient-patterns-2014-days180.txt", col_names = F,col_types = "cciii")
colnames(ref) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")

ref = ref[ref$NPI %in% short$NPI & ref$NPI2 %in% short$NPI,]

##Create sparse adjacency matrix
g = graph.edgelist(as.matrix(ref[,1:2]))
E(g)$weight = ref$Unique
W = get.adjacency(g, attr = 'weight', sparse = T)

# Regularize W
O = rowSums(W)
tau = rep(mean(O), length(O))

L_reg = Diagonal(length(O), 1/sqrt(O+tau))%*%W%*%Diagonal(length(O), 1/sqrt(O+tau))

gref = graph.adjacency(L_reg, weighted = T, diag = F)
core = coreness(gref)
hist(core)
gsm = induced_subgraph(gref,core>300)
M = get.adjacency(gsm,attr = "weight", sparse = T)

gref = graph.adjacency(L_reg, weighted = T, diag = F)
core = coreness(gref)
hist(core)
gsm = induced_subgraph(gref,core>median(core))
M = get.adjacency(gsm,attr = "weight", sparse = T)


##SpareAHC
hclust = sparseAHC(M, linkage = "average")

