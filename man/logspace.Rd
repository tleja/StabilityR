\name{logspace}
\alias{logspace}
\title{Logarithmically spaced values}
\description{Generate logarithmically spaced vector.}
\usage{
logspace(x1, x2, n) 
}
\arguments{
	\item{x1, x2}{The staring 10^x1 and (maximal) end 10^x2 values of the sequence.}
	\item{n}{Desired length of the sequence.}
}
\details{
The \code{logspace} function generates logarithmically spaced vector containing n values between 10^x1 and 10^x2.
}
\seealso{
	\code{\link{stability}}
}
\examples{
\dontrun{
# generate 50 values between 10^-3 and 10^1
logspace(x1=-3, x2=1, n=50)
}
}
