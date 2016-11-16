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
library(geosphere)

peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

short = read.csv("Data/NPI_zip.csv", colClasses = "character")

load(file = "Data/HClust_NPI.RData")

NPI = read.csv("Data/HClust_NPI.csv")


#cluster into HRR
k = 306
x = cutree(hclust,k)

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
print("Process started")
#Assigning zipcode to cluster
excl = zipcode[!(zipcode$zip %in% agg$zip), ]
excl$clus = 0
for(i in 1:dim(excl)[1]){
	tmp = excl[i,]
	loc = c(tmp$longitude, tmp$latitude)
    m = apply(cbind(agg$longitude, agg$latitude), 1, function(p, loc) distCosine(p,loc), loc = loc)
	excl$clus[i] = agg$clus[which.min(m)]
	if(i%%1000==0){print(i)}
}
write.csv(excl, "Data/AllZip_clus.csv", row.names = F)



