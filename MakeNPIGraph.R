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
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

short = read.csv("Data/NPI_zip.csv", colClasses = "character")
ref = read_csv(file = "Data/physician-shared-patient-patterns-2014-days180.txt", col_names = F,col_types = "cciii")
colnames(ref) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")

ref = ref[ref$NPI %in% short$NPI & ref$NPI2 %in% short$NPI,]

##Create sparse adjacency matrix
g = graph.edgelist(as.matrix(ref[,1:2]), directed = F)
E(g)$weight = ref$Unique
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


##SpareAHC
hclust = sparseAHC(M, linkage = "average")
save(hclust, file = "Data/HClust_NPI.RData")

NPI = as.data.frame(colnames(M))
colnames(NPI) = "NPI"
write.csv(NPI, "Data/HClust_NPI.csv", row.names = F)


NPI = read.csv("Data/HClust_NPI.csv")

load("Data/HClust_NPI.RData") #load hclust
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


# reorganize the clus
unique_label = sort(unique(agg$clus))
ind = c(1:length(unique(agg$clus)))
mapping = data.frame(unique_label,ind)

for(i in 1:length(mapping$unique_label)){
	agg$clus[agg$clus == mapping$unique_label[i]] = mapping$ind[i]
}

agg = merge(agg, zipcode, by = "zip", all.x = T)
loc =agg[,c("longitude", "latitude")]

##Mosaic plot
x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = agg$clus
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))

filename = paste("Pix/NPI_HRR", max(dt$clus), ".pdf", sep = "")
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





