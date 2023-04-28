\name{fieldboot}
\alias{fieldboot}

\title{Random Field Bootstrap}

\description{
The function bootstraps a data array representing a random field using the Moving Block Bootstrap, Circular Block Bootstrap or Nonoverlapping Block Bootstrap.
}

\usage{
	fieldboot(arr,func,B,blocklens,...,method=c("movingblock","nonoverlapping","circular"))
}

\arguments{
	\item{arr}{the data array}
	\item{func}{the function to apply to each sample.}
	\item{B}{the number of bootstrap samples.}
	\item{blocklens}{a scalar or a vector of block lengths.}
	\item{method}{method for array reconstruction.}
	\item{...}{optional additional arguments for the \code{func} function.}
}

\details{
The \code{arr} argument is a multidimensional array (with at least two dimensions).

The \code{blocklens} argument specifies the block lengths in all
dimensions. If it is a scalar integer value, the same size is used for
all dimensions, otherwise it must be an integer vector with as many elements
as there are dimensions in the array.

Different block bootstrap methods are supported to rebuild an array by bootstrapping
its blocks. The possible values of the \code{method} argument are: 
\emph{"movingblock"}, \emph{"nonoverlapping"}, or
\emph{"circular"}.
If it is not specified, the default method is \emph{"movingblock"}.
Method names may be abbreviated.

}

\value{
The \code{fieldboot} function returns an object of class \code{boodd} (see \link{class.boodd}).
}

\author{
	Bernard Desgraupes \cr
	\email{bernard.desgraupes@parisnanterre.fr}\cr
	University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{
P. Bertail, D. N. Politis, N. Rhomari (2000). Subsampling continuous parameter random fields and a Bernstein inequality,  \emph{Statistics}, 33(4), 367-392.

D.J. Nordman, S.N. Lahiri (2004). On optimal spatial subsample size for variance estimation, \emph{The Annals of Statistics}, 32(5), 1981-2027.

D.N.Politis, J.P.Romano (1993). Nonparametric Resampling for Homogeneous Strong Mixing Random Fields, \emph{J. Multivar. Anal.}, 47(2),   301-328.

D.N.Politis, J.P.Romano (1994). Large Sample Confidence Regions Based on Subsamples under Minimal Assumptions,  \emph{Ann. Statist.}, 22(4), 2031-2050.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.
}

\seealso{\code{\link{blockboot}},
\code{\link{jackVarBlock}},
\code{\link{jackVarField}}}

\examples{
# This example needs the RandomFields package
library(RandomFields)
# Simulate a gaussian random field with exponential covariance
arr <- model <- RMexp(var=1,scale=10)
# Simulate on [0,1]^2 a gid of size 100x100
from <- 0
to <- 10
x.seq <- seq(from, to, length=100)
y.seq <- seq(from, to, length=100)
arr <- RFsimulate(model,x=x.seq,y=y.seq)
arr <- as.array(arr)
# Bootstrap the array
B <- 299
blens <- c(20,10)
boo1 <- fieldboot(arr,mean,B,blens)
# Function with additional argument
meanpow <- function(arr,p) {mean(abs(arr)^p)^(1/p)}
boo2 <- fieldboot(arr,meanpow,B,blens,3)
# Circular blocks reconstruction
blens <- c(18,8)
boo3 <- fieldboot(arr,mean,B,blens,method="circular")
}

\keyword{random fields}