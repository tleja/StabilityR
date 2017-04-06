\name{graph_rings}
\alias{graph_rings}
\docType{data}
\title{Ring-of-rings graph}
\description{An undirected network containing 5 rings of 20 nodes in a ring configuration. The edges within the rings are 5 times stronger than the edges between the rings.}
\usage{
	data(graph_rings)
}
\format{
	\code{graph_rings} is a weighted adjacency matrix with 100 nodes.
}
\references{
	"Markov dynamics as a zooming lens for multiscale community detection: non clique-like communities and the field-of-view limit". Schaub, M. T., Delvenne, J.-C., Yaliraki, S. N. & Barahona, M. PLoS ONE, Public Library of Science, 2012, 7, e32210
}
\seealso{
	\code{\link{graph_demo}},
	\code{\link{forceAtlas2}},
	\code{\link{stability}},
}
\examples{
\dontrun{
	# load the graph
	data(graph_rings)

	# optimise the network layout with the forceAtas2 algorithm
	net = forceAtlas2(mat=graph_rings, iterations=4000, plotstep=10, k=100, ovd=10, seed=4)

	# visualise the network
	plotNetwork(mat=graph_rings, pos=net)
}
}
\keyword{datasets}
