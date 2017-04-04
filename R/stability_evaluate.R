stability_evaluate <- 
function(mat, res, cpu, opts) {
	
	suppressMessages(require(doSNOW, quietly=TRUE))
	
	if (is.null(res$call)) stop("res$call is missing")
	if (length(grep("evaluate",res$call))!=0) stop("evaluation was already performed")
	if (!is.matrix(mat)) stop("input should be a matrix") 
	if (!isSymmetric(mat)) stop("input matrix is asymmetric")
	if (nrow(mat)!=nrow(res$P)) stop("nrow(mat) != nrow(res$P)")
	if (!missing(cpu)) if (cpu < 1) stop("cpu must be an integer >= 1")
	
	# options to be passed to foreach
	if (missing(opts)) opts = NULL
	
	### prepare shared objects 
	par = res$call
	time = res$time
	p_unique = unique(res$P, MARGIN=2)
	if (min(p_unique)==0) p_unique = p_unique + 1
	
	PI = statDist(mat, type = par$type)
	
	### open parallel sockets
	nc = parallel::detectCores()
	if (!missing(cpu)) if (cpu < nc) nc = cpu
	cl = makeSOCKcluster(nc)
	registerDoSNOW(cl)
	# the clusterExport only to be used when sourced locally
	# parallel::clusterExport(cl, c("matExp",
	# 						'mat','par','time','p_unique','PI','res'))
	cat(paste('Number of registered CPU cores:', length(cl)), '\n')
	
	if (par$full) {
		cat('Evaluating stability: full,', par$type, '\n')
		tmp = foreach(i = seq_along(time), .options.snow=opts) %dopar% {
			solution = matExp(mat, time[i], type=par$type, prec=NULL)
			stabilities = evaluate_full(p=p_unique, s=solution, d=PI)
			ind = which.max(stabilities)
			s = stabilities[ind]
			p = p_unique[,ind]
			n = max(p)
			return(list(S=s, N=n, P=p))
		}
	}
	else {
		cat('Evaluating stability: linearised,', par$type, '\n')
		solution = mat/sum(mat)
		tmp = foreach(i = seq_along(time), .options.snow=opts) %dopar% {
			stabilities = evaluate_linearised(p=p_unique, s=solution, d=PI, time=time[i])
			ind = which.max(stabilities)
			s = stabilities[ind]
			p = p_unique[,ind]
			n = max(p)
			return(list(S=s, N=n, P=p))
		}
	}
	
	stopCluster(cl)
	
	res$S = do.call(c, lapply(tmp, function(x) x$S))
	res$N = do.call(c, lapply(tmp, function(x) x$N))
	res$P = do.call(cbind, lapply(tmp, function(x) x$P))
	
	res$call = call('stability_evaluate', mat=quote(mat), par)
	
	return(res)
}
