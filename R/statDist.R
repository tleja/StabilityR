statDist <-
function(mat, type='normalised') {
	if (!isSymmetric(mat)) stop("input matrix is asymmetric")
	type = match.arg(type, c('normalised','combinatorial'))
	mat = unname(mat)
	
	if (type=='normalised') {
		D = apply(mat,1, sum) # sum of weigths per node
		dist = D/sum(D)
	}
	if (type=='combinatorial') {
		N = ncol(mat) # node number
		dist = rep(1/N, N)
	}
	
	# todo: revise safety condition, sometimes it flags due to number rounding. 
	if (sum(dist)!=1) warning('sum(dist) != 1')
	return(dist)
}
