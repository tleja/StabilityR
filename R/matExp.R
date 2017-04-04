matExp <- 
function(mat, time, type = 'normalised', prec = 1e-9) {
	# require(rexpokit) # expokit_dgpadm_Qmat function
	
	if (!is.matrix(mat)) stop("input should be a matrix")
	if (!isSymmetric(mat)) stop("input matrix is asymmetric")
	# the adjacency matrix should be square (for directed, not implemented yet)
	# if (nrow(mat)!=ncol(mat)) stop("input matrix isn't square")
	type = match.arg(type, c('normalised','combinatorial'))
	
	colnames(mat) = rownames(mat) = NULL
	NbNodes = ncol(mat)
	sm = apply(mat,2,sum)
	
	if (type == 'normalised') {
		diagdeg = sm/sum(sm) # stationary distribution (normalised)
		trans = mat*(sm^(-1)) # (stochastic) transition matrix
		Lap = trans-diag(NbNodes)
		exponential = rexpokit::expokit_dgpadm_Qmat(Lap, t=time)
		solution = diagdeg*exponential
	}
	
	if (type == 'combinatorial') {
		av_degree = sum(mat)/NbNodes
		diagdeg = rep(1/NbNodes, NbNodes) # stationary distribution (combinatorial)
		Lap = -(mat-diag(sm))
		exponential = rexpokit::expokit_dgpadm_Qmat(Lap/av_degree, t=-time)
		solution = diagdeg*exponential
	}
	
	if (!is.null(prec)) {
		mx = apply(solution,2,max)
		solution = max(mx)*prec*round(solution/(max(mx)*prec))
	}
	
	return(solution)
}
