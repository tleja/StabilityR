.prim2 <- function(D) {
	
	# Number of nodes in the network
	N = nrow(D)
	
	# A matrix with the longest links in the paths of the MST
	LLink = matrix(0, ncol=N, nrow=N)
	
	# Allocate a matrix for the edge list 
	E = matrix(0, ncol=N, nrow=N)
	
	allidx = 1:N
	
	# Start with a node
	mstidx = 1
	otheridx = setdiff(allidx, mstidx)
	
	# A vector (T) with the shortest distance to all nodes.
	T = D[1,otheridx]
	P = rep(1, length(otheridx))
	
	while (length(T) > 0) {
		i = which.min(T)
		m = T[i]
		idx = otheridx[i]
		
		# Update the adjancency matrix	
		E[idx,P[i]] = 1
		E[P[i],idx] = 1
		
		# Update the longest links
		# 1) indexes of the nodes without the parent
		idxremove = which(mstidx == P[i])
		tempmstidx = mstidx
		tempmstidx = tempmstidx[-idxremove]
		
		# 2) update the link to the parent
		LLink[idx,P[i]] = D[idx,P[i]]
		LLink[P[i],idx] = D[idx,P[i]]
		
		# 3) find the maximal
		tempLLink = pmax(LLink[P[i],tempmstidx], D[idx,P[i]]);
		LLink[idx, tempmstidx] = tempLLink
		LLink[tempmstidx, idx] = tempLLink

		# As the node is added clear his entries		
		P = P[-i]	
		T = T[-i]
		
		# Add the node to the list 
		mstidx = c(mstidx , idx)

		# Remove the node from the list of the free nodes
		otheridx = otheridx[-i]

		# Updata the distance matrix
		Ttemp = D[idx, otheridx]
		
		if (length(T) > 0) {
			idxless = which(Ttemp < T)
			T[idxless] = Ttemp[idxless]
			P[idxless] = idx
		}
	}
	return(list(E=E, LLink=LLink))	
}

.rmst <- function(D, alpha) {
	N = nrow(D)

	mst = .prim2(D)
	
	# Find half the distance to nearest neighbors
	Dtemp = D + diag(N)*max(D)
	mD = apply(Dtemp,2,min)/alpha
	
	# Check condition
	mD = matrix(rep(mD,N),N, byrow=T) + matrix(rep(mD,N),N, byrow=F)
	E = (D - mD < mst$LLink)*1
	E = E - diag(diag(E))
	
	return(E)
}

.addToNetwork <- function(Danchor, Dnew, alpha) {
	if (isSymmetric(Danchor)!=TRUE) stop('data matrix is not symmetric: distance matrix is required')
	
	Eanchor = .rmst(Danchor, alpha)
	LLinkanchor = .prim2(Danchor)$LLink
	
	N = nrow(Eanchor)
	Nnew = nrow(Dnew)	
	
	if (N != ncol(Dnew)) stop('ncol(Dnew) must match to dim(Danchor)')
	
	# Find the nearest neighbor of each anchor point
	Danchor = Danchor + diag(N)*max(Danchor) 
	nnanchor = apply(Danchor,2,min)
	
	# Edges between the anchor points and the new points
	Enew = matrix(0, nrow=N, ncol=Nnew)
	
	for (i in 1:Nnew) {
		Dp = Dnew[i,]

		# Find the nearest neighbor
		nn = which.min(Dp)
		m = Dp[nn]
		# message('Nearest Neighbor : ', nn);	
		
		# Add an edge between the new point and the nearest neighbor
		Enew[nn, i] = 1	

		# Find the longest link between the nearest neighbor and the other anchor points
		tempLLink = LLinkanchor[nn,]

		tempLLink = pmax(tempLLink, m);

		idx = which(Dp - (m+nnanchor)/alpha < tempLLink);
		# message('Connecting to ',length(idx),' more')	
		Enew[idx, i] = 1
	}
	
	return(Enew)	
}

.getDist <- function(mat) {
	# correlation matrix
	if (all(diag(mat)==1)) {
		if (max(mat)>1 | min(mat)<(-1)) stop('matrix is not a correlation matrix')
		
		# correlation matrix to distance (0-1), where (1 => 0) and (-1 or floor(min) => 1)
		mmin = floor(min(mat))
		d = (1 - mat) / (1 - mmin)	
	}
	# distance matrix
	if (all(diag(mat)==0)) {
		if (min(mat)<0) stop('matrix is not a distance matrix')
		
		# set maxiumum range 0-1 for distance matrix
		mmax = ceiling(max(mat))
		d = mat/mmax
	}
	return(d)
}


### MAIN function
RMST <-
function(mat, grow, alpha=2, allEdges=TRUE, verbose=FALSE) {
	if (!is.matrix(mat)) stop("input should be a matrix")
	if (!isSymmetric(mat)) stop("input matrix is asymmetric")
	
	if (missing(grow)) {
		d = .getDist(mat)
		# apply rmst and multiply by the inverted distance; where 1 = highest weight
		s = .rmst(d, alpha=alpha)
		sparse = s*(1-d)
	}
	else {
		if (is.null(rownames(mat)) | is.null(colnames(mat))) stop('mat rownames or colnames are empty')
		if (!identical(rownames(mat), colnames(mat))) stop('mat rownames != colnames')
		if (!is.list(grow)) stop('grow must be provided in a list format')
		if (length(grow)>2) stop('only grow of length 2 is allowed')
		
		g1 = intersect(sort(unique(grow[[1]])), rownames(mat))
		g2 = intersect(sort(unique(grow[[2]])), rownames(mat))
		if (length(intersect(g1,g2))>0) stop('grow contains overlapping nodes')
		subs = c(g1, g2)
		
		if (verbose) cat(sprintf('N1: %s \nN2: %s \nN12: %s', length(g1), length(g2), length(subs)))
		
		# get distance matrices for g1, g2 and merged
		d = .getDist(mat[subs,subs])
		d1 = d[g1,g1]
		dm = d[g2,g1] # g2 in rows and g1 in columns
		if (allEdges) d2 = d[g2,g2]
		
		# apply rmst and addToNetwork
		s1 = .rmst(d1, alpha=alpha)
		sm = .addToNetwork(d1, dm, alpha=alpha)
		if (allEdges) s2 = .rmst(d2, alpha=alpha)
	
		# save the results to a matrix
		m = matrix(0, nrow=length(subs), ncol=length(subs))
		rownames(m) = colnames(m) = subs
		
		m[g1,g1] = s1
		m[g1,g2] = sm
		m[g2,g1] = t(sm)
		if (allEdges) m[g2,g2] = s2
		if (!isSymmetric(m)) warning('rmst matrix is asymetric')
		
		sparse = m*(1-d)
		sparse = sparse[order(subs),order(subs)]
		
		type = rep(NA, nrow(sparse))
		type[rownames(sparse) %in% g1] = 1
		type[rownames(sparse) %in% g2] = 2
		
		sparse = list(mat=sparse, type=type)
	}
	
	return(sparse)
}
