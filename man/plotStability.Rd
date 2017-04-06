\name{plotStability}
\alias{plotStability}
\alias{plotStability_VI}
\alias{plotStability_heatmap}
\title{Stability plotting functions}
\description{Plotting HTML functions for \code{stability} output.}
\usage{
plotStability(res)
plotStability_VI(res)
plotStability_heatmap(res)
}
\arguments{
	\item{res}{The output from the \code{\link{stability}} function.}
}
\details{
The \code{plotStability} function visualises the stability values with the corresponding number of partitions at each Markov time. 

The \code{plotStability_VI} shows the variation of information values obtained from M iterations of stability optimisation at each Markov time. It shows the robustness of partitions at each Markov time.

The \code{plotStability_heatmap} displays pair-wise variation of information calculated between the partition vectors obtained across the Markov time scale. It shows the similarities between the resulting partitions.
}
\seealso{
	\code{\link{stability}}
}
\examples{
\dontrun{
# load the graph
data(graph_rings)

# execute the stability function
mat = graph_rings
time = logspace(-3,1,61)
res1 = stability(mat=mat, time=time, full=TRUE, type='norm', M=100)
res2 = stability_evaluate(mat=mat, res=res1)

# visualise the stability values and number of partitions
plotStability(res2)

# visualise the robustness of partitions at each Markov time
plotStability_VI(res2)

# compare the resulting partitions across the Markov time scale
plotStability_heatmap(res2)
}
}
