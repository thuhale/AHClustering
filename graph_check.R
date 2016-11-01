library(Matrix)
library(sparseAHC)
library(igraph)
library(vegan)

A = cbind(c(0,12,0,0,0,0,14,0,0),
		  c(12,0,1,11,0,0,13,0,0),
		  c(0,1,0,3,0,0,0,0,0),
		  c(0,11,3,0,2,5,0,0,0),
		  c(0,0,0,2,0,10,0,9,0),
		  c(0,0,0,5,10,0,6,7,0),
		  c(14,13,0,0,0,6,0,0,4),
		  c(0,0,0,0,9,7,0,0,8),
		  c(0,0,0,0,0,0,4,8,0))
rownames(A) = c("A", "B", "C", "D", "E", "F", "G", "H","I") 
colnames(A) = rownames(A)

#par(mfrow = c(1,2))
g = graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
E(g)$label = E(g)$weight
E(g)$id <- seq_len(ecount(g))
mst = mst(g, weight = 1/E(g)$weight)
E(g)$color <- "grey"
E(g)$color[E(mst)$id] <- "yellow"
plot(g, edge.width = E(g)$weight)


D = rowSums(A)
tau = rep(mean(D), length(D)) 
L = diag(1/sqrt(D+tau))%*%A%*%diag(1/sqrt(D+tau))
colnames(L) = colnames(A)
g = graph_from_adjacency_matrix(L, mode = "undirected", weighted = TRUE)
E(g)$label = round(E(g)$weight,4)*100
E(g)$id <- seq_len(ecount(g))
mst = mst(g, weight = 1/E(g)$weight)
E(g)$color <- "grey"
E(g)$color[E(mst)$id] <- "yellow"
plot(g, edge.width = E(g)$weight*20)

par(mfrow = c(1,2))
plot(sparseAHC(as(A,"sparseMatrix"), linkage = "single"), labels = colnames(A))
plot(sparseAHC(as(L,"sparseMatrix"), linkage = "single"), labels = colnames(L))






