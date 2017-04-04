edgeList <-
function(mat, sparse=TRUE, diag.rm=FALSE, verbose=TRUE) {
	# require(reshape2)
	
	if (!is.matrix(mat)) stop("input should be a matrix")
	if (nrow(mat)!=ncol(mat)) stop("input matrix isn't square")
	if (!identical(rownames(mat),colnames(mat)))
		stop("row names do not match with column names in the input matrix")
	
	if (isSymmetric(mat)) {
		if (verbose) cat("Processing undirected graph","\n")
		if (all(diag(mat)==0) | diag.rm) {
			if (verbose) cat("... removing matrix diagonal","\n")
			mat[upper.tri(mat, diag=T)] = NA
		}
		else mat[upper.tri(mat, diag=F)] = NA
		d = na.omit(reshape2::melt(mat, varnames=c("row","col"), value.name='weight'))
	}
	else {
		if (verbose) cat("Processing directed graph","\n")
		if (all(diag(mat)==0) | diag.rm) {
			if (verbose) cat("... removing matrix diagonal","\n")
			diag(mat) = NA
		}
		d = na.omit(reshape2::melt(t(mat), varnames=c("col","row"), value.name='weight'))
		d = data.frame(from=d$row, to=d$col, weight=d$weight, stringsAsFactors=FALSE)
	}
	
	if (sparse) {
		if (verbose) cat("... removing null edges","\n")
		d = d[d$weight!=0,]
		rownames(d) = NULL
	}
	
	return(d)
}
