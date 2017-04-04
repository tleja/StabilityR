stability <- 
function(mat, time, full = FALSE, type = 'normalised',
					M = 100, VI = TRUE, seed = NULL, cpu, prec = 1e-9, opts, log = FALSE) {
	
	# load parallel packages
	suppressMessages(require(doSNOW, quietly=TRUE))
	
	### evaluate laplacian type to be computed
	type = match.arg(type, c('normalised','combinatorial'))
	
	### evaluate input adjacency matrix and time vector
	if (!is.matrix(mat)) stop("input should be a matrix") 
	# if (nrow(mat)!=ncol(mat)) stop("input matrix isn't square")
	if (!isSymmetric(mat)) stop("input matrix is asymmetric")
	if (!identical(rownames(mat),colnames(mat)))
		stop("row names do not match with column names in the input matrix")
	if (!is.vector(time) | !is.numeric(time)) stop("time must be a numeric vector")
	if (!missing(cpu)) if (cpu < 1) stop("cpu must be an integer >= 1")
	
	### options for foreach
	if (missing(opts)) opts = NULL
	
	### save node names	
	nodes = colnames(mat)
	rownames(mat) = colnames(mat) = NULL
	
	### initiate timestamp
	strt = Sys.time()
	
	### compute stationary distribution	
	cat(paste('Computing stationary distribution:', type, '...'))
	d = statDist(mat, type=type)
	dist = cbind(d,d)
	cat(paste(' done in', round(as.numeric(Sys.time()-strt, 
							unit='mins'), 3), 'mins'), '\n')
	
	### create log for debugging
	logfile = paste0("log_stability_", format(Sys.time(), "%Y%m%d_%H%M%S"))
	if (log) writeLines(c(""), logfile)

	### record number of CPU cores to be used
	nc = parallel::detectCores()
	if (!missing(cpu)) if (cpu < nc) nc = cpu
	cl = makeSOCKcluster(nc)
	registerDoSNOW(cl)
	# the clusterExport only to be used when sourced locally
	# parallel::clusterExport(cl, c("matExp","edgeList","run_gen_louvain2",
	# 					"readStability","VI","varinfo",
	# 					'mat','time','dist','log'))
	cat(paste('Number of registered CPU cores:', length(cl)), '\n')
	
	### execute stability (parallelised)
	if (full) {
		cat('Computing stability: full ...')
		res = foreach(i = seq_along(time), .options.snow=opts) %dopar% {
			
			if (log) cat(paste("... stability at t[ind] =",i,'at', 
						round(as.numeric(Sys.time()-strt, unit='mins'), 3), 'mins','\n'), 
							file=logfile, append=TRUE)
			
			### compute matrix exponential and input graph
			m = matExp(mat, time=time[i], type=type, prec=prec)
			if (!isTRUE(all.equal(sum(m),1))) 
				warning(paste("sum of matrix exponential != 1 for t[ind] =", i))
			e = edgeList(m, sparse=T, verbose=F)
			graph = cbind(e$col-1, e$row-1, e$weight)
			
			### run_gen_louvain	and replicate over M	
			if (!is.null(seed)) set.seed(seed)	
			out = run_gen_louvain(graph=graph, nullmod=dist, time=1, M=M)

			### output processing
			out = readStability(out, VI=VI)
			
			### cleaning files and output
			rm(graph)
			return(out)
		}
		cat(paste(' done in', round(as.numeric(Sys.time()-strt, 
								unit='mins'), 3), 'mins'), '\n')
	}
	else {
		cat('Computing stability: linearised ...')

		### compute input graph
		e = edgeList(mat/sum(mat), sparse=T, verbose=F)
		graph = cbind(e$col-1, e$row-1, e$weight)
		
		### run_gen_louvain
		res = foreach(i = seq_along(time), .options.snow=opts) %dopar% {

			if (log) cat(paste("... stability at t[ind] =",i,'at', 
						round(as.numeric(Sys.time()-strt, unit='mins'), 3), 'mins','\n'), 
							file=logfile, append=TRUE)

			### run_gen_louvain and replicate over M
			if (!is.null(seed)) set.seed(seed)	
			out = run_gen_louvain(graph=graph, nullmod=dist, time=time[i], M=M)

			### output processing
			out = readStability(out, VI=VI)

			### file output
			return(out)
		}
		cat(paste(' done in', round(as.numeric(Sys.time()-strt, 
								unit='mins'), 3), 'mins'), '\n')
	}
	
	stopCluster(cl)	

	### compute final statistics
	cat('Computing final statistic ...')
	S = do.call(c, lapply(res, function(x) x$S))
	N = do.call(c, lapply(res, function(x) x$N))
	V = do.call(c, lapply(res, function(x) x$VI))
	P = do.call(cbind, lapply(res, function(x) x$P))
	cat(paste(' done in', round(as.numeric(Sys.time()-strt, 
							unit='mins'), 3), 'mins'), '\n')	
	
	### save call parameters
	par = call('stability', mat=quote(mat), 
				time=substitute(logspace(x1,x2,n=n), 
						list(x1=log10(min(time)), x2=log10(max(time)), 	
								n=as.numeric(length(time)))), 
				full=full, type=type, M=M, VI=VI, seed=seed, cpu=as.numeric(nc))
	
	### return stability results in a list format
	return(list(S=S, N=N, VI=V, P=P, names=nodes, time=time, call=par))
}
