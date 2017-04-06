\name{stability_evaluate}
\alias{stability_evaluate}
\title{Stability evaluation}
\description{The \code{stability_evaluate} function evaluates stability function against the observed graph partitions.}
\usage{
stability_evaluate(mat, res, cpu, opts) 
}
\arguments{
	\item{mat}{Adjacency matrix of an undirected graph. See \code{\link{stability}}}
	\item{res}{List containing the ouput from the \code{\link{stability}} function.}
	\item{cpu}{Number of CPU nodes to be used for parallel execution. The function automatically detects the number of available nodes, while \code{cpu=1} will execute it in sequential mode.}
	\item{opts}{Optional list containing functions to be passed during parallel execution. See \code{\link{stability}}}
}
\details{
The \code{stability_evaluate} function uses the output from the \code{stability} function to refine the stability values and the optimal partitions.
}
\value{
See \code{\link{stability}}.
}
\seealso{
	\code{\link{stability}},
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
res1 = stability(mat=mat, time=time, full=TRUE, type='norm', M=100)

# perform stability evaluation
res2 = stability_evaluate(mat=mat, res=res1)

# visualise the results
plotStability_heatmap(res1)
plotStability_heatmap(res2)
}
}
