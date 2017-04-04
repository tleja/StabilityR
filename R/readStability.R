readStability <- 
function(x, VI=TRUE) {
	p = x$partition
	s = x$stability
	v = 0
	if (VI) v = VI(p, out='mean')
	p = p[,which.max(s)] + 1
	n = max(p)
	s = max(s)
	return(list(S=s, N=n, VI=v, P=p))
}
