\name{stability}
\alias{stability}
\title{Stability of a graph partition}
\description{The \code{stability} function performs community detection using the stability of a graph partition.}
\usage{
stability(mat, time, full = FALSE, type = "normalised", M = 100, 
    VI = TRUE, seed = NULL, cpu, opts, log = FALSE, prec = 1e-09) 
}
\arguments{
	\item{mat}{
	Adjacency matrix of an undirected graph.
	}
	\item{time}{
	Numeric vector with stopping, Markov time values. See \code{\link{logspace}}.
	}
	\item{full}{
	If \code{TRUE}, the function performs calculation of the full stability, while \code{FALSE} executes linearised stability (default).
	}
	\item{type}{
	Type of the laplacian that should be used to calculate stability. Options include "normalised" (default) or "combinatorial".
	}
	\item{M}{
	Number of optimisations to be performed at each Markov time. Used to improve robustness of the partitions and to calculate the variation of information.
	}
	\item{VI}{
	If \code{TRUE}, the variation of information is calculated (default).
	}
	\item{seed}{
	Optional value to define starting seed for the louvain optimisations at each Markov time. Used to obtain reproducible results.
	}
	\item{cpu}{
	Number of CPU nodes to be used for parallel execution. The function automatically detects the number of available nodes, while \code{cpu=1} will execute it in sequential mode.
	}
	\item{opts}{
	Optional list containing functions to be passed during parallel execution. See examples below.
	}
	\item{log}{
	If \code{TRUE} the ouput log file is created which allows to track the process of the stability function.
	}
	\item{prec}{
	Precision argument, defines a threshold for the range of weights allowed in the laplacian exponential matrix of the full stability.
	}
}
\details{
The \code{stability} function finds the optimal partitions of the undirected graph by optimising the stability at each Markov time.
}
\value{
The function returns a list with the following elements:
	\item{S}{Vector with stability values at each Markov time.}
	\item{N}{Vector with number of partitions detected at each Markov time.}
	\item{VI}{Vector with variation of information values calculated at each Markov time.}
	\item{P}{Matrix combining partition vectors obtained at each Markov time.}
	\item{names}{Character vector with names of nodes extracted from the input graph. \code{NULL} if unnamed graph was provided.}
	\item{time}{Vector with Markov time values used for stability optimisation.}
	\item{call}{Object of mode \code{\link{call}} containing an unevaluated stability function call with arguments used during the execution.}
}
\references{
"Stability of graph communities across time scales", Delvenne, J.-C., Yaliraki, S. N. & Barahona, M. arXiv:0812.1811 (2009) and then published in Proceedings of the National Academy of Sciences, 2010, 107, 12755-12760;

"Laplacian Dynamics and Multiscale Modular Structure in Networks" Lambiotte, R., Delvenne, J.-C. & Barahona, M. arxiv:0812.1770 (2009)

To optimize the stability quality function, we use the Louvain algorithm as described: "Fast unfolding of communities in large networks", Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre, Journal of Statistical Mechanics: Theory and Experiment 2008 (10), P10008
}
\seealso{
	\code{\link{stability_evaluate}},
	\code{\link{plotStability}}
}
\examples{
\dontrun{
# load the graph
data(graph_rings)
mat = graph_rings

# set the Markov time range
time = logspace(-3,1,61)

# execute the full stability using normalised laplacian
res = stability(mat=mat, time=time, full=TRUE, type='norm', M=100)

# execute the function in sequential mode
res = stability(mat=mat, time=time, full=TRUE, type='norm', M=100, cpu=1)

# set up a progress bar using opts argument
pb = txtProgressBar(max=length(time), style=3)
opts = list(progress = function(n) setTxtProgressBar(pb, n))
res = stability(mat=mat, time=time, full=TRUE, type='norm', M=100, opts=opts)
close(pb)

# visualise the results
plotStability(res)
plotStability_VI(res)
plotStability_heatmap(res)
}
}
