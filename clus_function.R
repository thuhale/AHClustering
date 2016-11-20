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
