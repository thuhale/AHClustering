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
plot(sparseAHC(as(A,"sparseMatrix"), linkage = "single"), labels = colnames(A), main = "AHC on A")
plot(sparseAHC(as(L,"sparseMatrix"), linkage = "single"), labels = colnames(L), main = "AHC on L")


#Search the tree


slideTree = function(g){
	tau = 1/length(V(g))
	best = 1000
	to_cut = 1
	E(g)$cut = rep(best, length(E(g)))
	for(i in 1:length(E(g))){
		g1 = g
		g1 = delete_edges(g1, E(g1)[i])
		
		comp = components(g1)$csize
		cut = E(g)$weight[i]/((comp[1]+tau)*(comp[2]+tau))
				
		if(cut<best){
			best = cut
			to_cut = i
		}
		E(g)$cut[i] = cut
	}
	g1 = g
	g1 = delete_edges(g1, E(g1)[to_cut])

	mem = components(g1)$membership
	left_tree = induced_subgraph(g, V(g)[mem==1])
	right_tree = induced_subgraph(g, V(g)[mem==2])
	print(V(left_tree))
	print(V(right_tree))
	return(list(left_tree, right_tree))
}


x = slideTree(mst)
left_tree = x[[1]]
right_tree = x[[2]]

cutTree = function(g, dt){
	if(length(V(g)) == 1){
		return(g)
	}
	x = slideTree(g)
	left_tree = x[[1]]
	right_tree = x[[2]]
	cutTree(left_tree)
	cutTree(right_tree)
	
}



cutTree(mst)

## reorder dendogram









