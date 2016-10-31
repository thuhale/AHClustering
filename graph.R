rm(list = ls())
library(igraph)
library(readr)
library(Matrix)
library(sparseAHC)
library(MCL)
library(rARPACK)

library(zipcode)
data(zipcode) 
library(data.table)
zipcode = as.data.table(zipcode); 
setkey(zipcode,zip)
library(maps)
library(RColorBrewer)
library(randomcoloR)

compare = read.table("compare.csv", sep = ",", header=T, colClasses = "character")
PCP = compare$NPI[compare$PCP==1]
df = read.table("Data/PCP.csv", sep = ",", header=T, colClasses = "character")
d = df[!(df$zip %in% zipcode$zip),]


refer2011 = read_csv(file = "Data/physician-shared-patient-patterns-2011-days30.txt", col_names = F,col_types = "cciii")

colnames(refer2011) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")
refer2011 = refer2011[refer2011$NPI %in% compare$NPI & refer2011$NPI2 %in% compare$NPI,]


g = graph.edgelist(as.matrix(refer2011[,1:2]))
E(g)$weight = refer2011$Ties

#SPARSE AHC
W = get.adjacency(g, attr = 'weight', sparse = T)
O = rowSums(W)
tau = rep(1, length(O))

L_reg = Diagonal(length(O), 1/sqrt(O+tau))%*%W%*%Diagonal(length(O), 1/sqrt(O+tau))

L = L_reg[rownames(L_reg) %in% PCP,]
L = L[, !(colnames(L_reg) %in% PCP)] ##L is the bi-partite matrix
L = L%*%t(L) 
diag(L) = 0 ## remove self-loop

gref = graph.adjacency(L, weighted = T, diag = F)
core = coreness(gref)
hist(core)

gsm = induced_subgraph(gref,core>100)
sum(degree(gsm))
A = get.adjacency(gsm,attr = "weight")

diag(A) = 0
dim(A)
hc = sparseAHC(A, linkage = "average")
k = 400
clus = cutree(hc,k)

tab = as.data.frame(table(clus))
colnames(tab) = c("clus", "size")
tab = tab[sort(tab$size, decreasing = T, index.return = T)$ix,]
big = tab[tab$size>100,]
big_clus = big$clus
big_clus = as.numeric(paste(big_clus))

pal = distinctColorPalette(length(big_clus))
pal = data.frame(pal, big_clus)
colnames(pal) = c("col", "clus")

dt = data.frame(colnames(A), clus)
colnames(dt) = c("NPI", "clus")
dt = dt[dt$clus %in% big_clus,]
dt = merge(dt, df, by = "NPI", all.x = T)
dt = merge(dt, pal, by = "clus", all.x = T)

loc = zipcode[dt$zip,c("longitude", "latitude"), with= F]
map("state", col = 1)
points(loc, col = dt$clus, pch = 19, cex = 0.1)
i = 5
points(loc[dt$clus==big_clus[i]], col = dt$col[dt$clus==big_clus[i]], pch = 19, cex = 0.1)
title(main="SparseAHC")


#cluster
order = hc$order
NPI = colnames(A)[order[1:100]]
doc = data.frame(NPI, c(1:100))
colnames(doc) = c("NPI", "rank")

doc = merge(doc, df, by = "NPI", all.x = T)
doc = doc[sort(doc$rank, index.return = T)$ix,]

loc = zipcode[doc$zip,c("longitude", "latitude"), with= F]
map("state", col = 1)
points(loc[1:100], col = "blue", pch = 19, cex = 0.1)
title(main="SparseAHC")

city = zipcode[doc$zip, c("city", "state"), with = F]

