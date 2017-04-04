VI <- 
function(p.mat, out=c('mean','all')) {
	out = match.arg(out)
	
	vim = varinfo(p.mat+1)
	
	if (out=='mean') return(mean(vim[lower.tri(vim)]))
	if (out=='all') return(vim)
}
