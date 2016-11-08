rm(list = ls())
library(igraph)
setwd("/Users/lehathu/Desktop/AHC")
library(Matrix)
library(sparseAHC)
library(MCL)
library(zipcode)
data(zipcode)
library(maps)
library(deldir)
library(RColorBrewer)
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-US continential 
zipcode = zipcode[!(zipcode$state %in% peri),]
## This file reads in PCP edge_list network "PCP_Referral_List.csv" generated from MakePCPGraph.R file and performing cluster on that. 
## It also plots the cluster on map. Doctor's location is found on "NPI_zip.csv"
# read in doctor's zipcode
short = read.csv("Data/NPI_zip.csv", colClasses = "character")
#create a merge matrix
PCPRef = read.csv("Data/PCP_Referral_List.csv")
PCPRef$NPI = as.character(PCPRef$NPI)
PCPRef$NPI2 = as.character(PCPRef$NPI2)
g = graph.edgelist(as.matrix(PCPRef[,1:2]))
E(g)$weight = PCPRef$Weight
M = get.adjacency(g, attr = 'weight', sparse = T)

hclust = sparseAHC(M, linkage = "average")
order = order.dendrogram(as.dendrogram(hclust))
order = sort(order, index.return = T)$ix

sorted_order = data.frame(colnames(M), 1:length(colnames(M)), order)
colnames(sorted_order) = c("NPI", "index", "order")
sorted_order = merge(sorted_order, short, by = "NPI", all.x = T)
sorted_order = merge(sorted_order, zipcode, by = "zip", all.x = T)
sorted_order = sorted_order[,1:6]
sorted_order = sorted_order[order(sorted_order$index),]

merge = hclust$merge
merge = as.data.frame(merge)
colnames(merge) = c("V1", "V2")

write.csv(sorted_order, "Data/sorted_order.csv", row.names = F)
write.csv(merge, "Data/merge.csv", row.names = F)

#Point-wise plot
k = 500
x = cutree(hclust, k)
dt = data.frame(colnames(M), x)
colnames(dt) = c("NPI", "clus")
dt = merge(dt, short, by = "NPI", all.x = T)
dt = merge(dt, zipcode, by = "zip", all.x = T)


map("state", col = 1)
for(i in 1:k){
	points(dt$longitude[dt$clus==i], dt$latitude[dt$clus==i], col = i, pch = 19, cex = 0.5)
}

##Mosaic plot

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
#tab = as.data.frame(table(agg$clus))
#colnames(tab) = c("clus", "size")
#tab1 = tab[tab$size>40,]
#agg = agg[agg$clus %in% tab1$clus,]

# find a coloring of the graph g so that adjacent nodes have a different color.
# think of each "region" as a node, we want to draw a map of regions such that boundaries
#  are easily visible.
color.igraph = function(g){
 n = length(V(g))
 ordrank = order(-degree(g))
 neigh  = neighborhood(graph = g,order = 1,mode = "all")
 possibleCols = matrix(T,nrow = n, ncol = max(degree(g)))
 col = rep(NA,n)
 for(i in 1:n){
   id = ordrank[i]
   # find color
   col[id] = which(possibleCols[id,])[1]
   # make adjacent not have this color
   possibleCols[neigh[[id]],col[id]] = F
 }  

 brew = brewer.pal(max(col), "Paired") 
 return(brew[col])
}

aggregate.igraph = function(g, cw, edgelist =F){
 # g is a graph.
 # cw is a crosswalk, an nx2 data.table; first column is node id in g.  second column is new node id. the key should be the first column.
 # the number of unique entries in column 2 should be << n.  
 # if edgelist = T, then g is an edgelist.

 # This code doesn't work..
 # cw = as.data.table(cw)
 #   setkey(cw,colnames(cw)[1])
 #   

 if(!is.vector(cw)){
   setnames(cw,1:2, c("V1","V2"))
   setkey(cw,"V1")

   if(edgelist) el = g
   if(!edgelist) el = get.edgelist(graph = g)

   Enew = cbind( cw[el[,1]]$V2,
                 cw[el[,2]]$V2
   )  
   #   Enew = cbind( cw[as.character(as.numeric(el[,1]))]$V2,
   #                 cw[as.character(as.numeric(el[,2]))]$V2
   #   )
 }
 if(is.vector(cw)){
   if(edgelist) el = g
   if(!edgelist) el = get.edgelist(graph = g)

   Enew = cbind( cw[el[,1]],
                 cw[el[,2]]
   )  
 }
 if(edgelist) {if(ncol(el)>2) Enew = cbind(Enew, el[,3])}
 if(!edgelist) {if(length(E(g)$weight)>0) Enew = cbind(Enew, E(g)$weight)}
 if(ncol(Enew) != 3)  Enew = cbind(Enew, rep(1, nrow(Enew)))


 colnames(Enew)[3] = "weight"
 Enew = Enew[complete.cases(Enew),]
 gnew =graph.edgelist(el = Enew[,1:2],directed = T)
 E(gnew)$weight = as.numeric(Enew[,3])

 return(simplify(gnew, remove.loops= F))
}




# pass a dl object and a coloring
plot.dl = function(dl, clust, res = 1){
 # dl is a deldir tesselation
 # clust is a cluster assignment
 # res is resolution in degrees (could be improved) 
 # 1 degree is 70 miles. 

 #make a graph of adjacent clusters
 # make graph of zips
 # aggregate with clust
 el = dl$dirsgs[,5:6]
 g = graph.edgelist(as.matrix(el), directed =  F)
 gclust = simplify(aggregate.igraph(g,clust))
 # 

 # make a coloring
 cols = color.igraph(gclust)

 # find lines which separate clusters
 seps = apply(dl$dirsgs[,5:6],1, function(x, clust) return(clust[x[1]] != clust[x[2]]), clust)
 sdir = dl$dirsgs[seps,]
 # find lines which are less than res
 lengthLines = sqrt((sdir[,1] - sdir[,3])^2 + (sdir[,2] - sdir[,4])^2)
 map("state")

 for(clusterLabel in 1:max(clust-1)){
 good= (clust[sdir[,5]] == clusterLabel | clust[sdir[,6]] == clusterLabel)  & !sdir[,7] & ! sdir[,8]
 pv = as.matrix(sdir[good,1:4])
 el = t(apply(pv,1, function(x) return(c(paste(x[1],x[2]), paste(x[3],x[4])))))
 hashVec = rnorm(2)

 vertexLoc = rbind(pv[,1:2],pv[,3:4])
 nodehash = vertexLoc%*%hashVec
 uh = unique(nodehash)
 numel = cbind(match(pv[,1:2]%*%hashVec, uh),match(pv[,3:4]%*%hashVec,uh))
 vertexLoc = vertexLoc[match(uh,nodehash),]

 gpoly = graph.edgelist(numel,directed = F)
 # V(gpoly)$latlon = V(gpoly)$name
 # V(gpoly)$name=1:length(V(gpoly))
 cc = clusters(gpoly)
 polyList = list()
 polyList[[1]] = cc$no
 for(i in 1:cc$no){
   starting = which(cc$membership==i)[1]
   path = graph.dfs(gpoly,starting, unreachable = F)$order
   polyList[[i+1]] = path[complete.cases(path)]
 }

 for(pp in 1: polyList[[1]]){
   tmp = vertexLoc[polyList[[pp+1]],]
   polygon(tmp,col = cols[clusterLabel])  
 }
 }
 map("state", add=T)
}


agg = merge(agg, zipcode, by = "zip", all.x = T)
loc =agg[,c("longitude", "latitude")]

x = rnorm(loc$longitude,loc$longitude,.00001); y = rnorm(loc$latitude,loc$latitude, .00001)
clust = agg$clus
# dl = deldir(x,y)
dl = deldir(c(x,-1000,-1000,1000,1000),c(y,-1000,1000,1000,-1000))
#clust = c(clust, rep(max(clust)+1,4))


pdf(file = "Regions.pdf", height = 50, width = 100)
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



