\name{graph_demo}
\alias{graph_demo}
\docType{data}
\title{Multiscale graph}
\description{A simple undirected network with a multiscale community structure.}
\usage{
	data(graph_demo)
}
\format{
	\code{graph_demo} is a weighted adjacency matrix with 640 nodes.
}
\seealso{
	\code{\link{graph_rings}},
	\code{\link{forceAtlas2}},
	\code{\link{stability}}
}
\examples{
\dontrun{
	# load the graph
	data(graph_demo)

	# optimise the network layout with the forceAtas2 algorithm
	# it may take a few minutes for a large network
	net = forceAtlas2(mat=graph_demo, iterations=1000, plotstep=10, k=100, gravity=1, seed=4, ovd=10)

	# visualise the network
	plotNetwork(mat=graph_demo, pos=net)
}
}
\keyword{datasets}
