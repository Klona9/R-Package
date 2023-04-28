\name{sieveboot}
\alias{sieveboot}

\title{Autoregressive Sieve Bootstrap}

\description{
The function applies autoregressive sieve bootstrap to stationary time
series. The idea is to estimate an \emph{AR(p)} model, with p large and to
resample the centered estimated residuals to reconstruct an \emph{AR(p)}
bootstrap times series on which a given statistics is computed.
}

\usage{
	sieveboot(x,func,B,order=NULL,...)
}

\arguments{
	\item{x}{a vector of observations.}
	\item{func}{the function to apply to each sample.}
	\item{B}{a positive integer: the number of bootstrap samples.}
	\item{order}{a positive integer: the order of the autoregressive process.}
	\item{...}{optional additional arguments for the \code{func} function.}
}

\details{
The sieve bootstrap procedure approximates the given data by a large order autoregressive process.

The validity of the sieve bootstrap was obtained for \eqn{p=o((n/log(n))^(1/4))}. 
When the order is not specified, we set
\eqn{p=4*(n^1/4)/log(n)^(1/2)}.
}

\value{
\code{sieveboot} returns an object of class \code{boodd} (see \link{class.boodd}).
}

\author{
	Bernard Desgraupes \cr
	\email{bernard.desgraupes@parisnanterre.fr}\cr
	University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\source{Bühlmann, P. (1997). Sieve Bootstrap for time series. Bernoulli, 3, 123-148.}

\references{
Bühlmann, P. (1997). Sieve Bootstrap for time series.
\emph{Bernoulli,} \bold{3}, 123-148.

Choi, E., Hall, P. (2000). Bootstrap confidence regions computed
from autoregression of abitrary order, \emph{Journal of the Statistical
Royal Society}, ser. B, \bold{62}, 461-477.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{freqboot}},
\code{\link{aidedboot}}}

\examples{
n <- 200
B <- 299
x <- arima.sim(list(order=c(0,0,4),ma=c(0.7,0.4,-0.3,-0.1)),n=n)
sieveboot(x,mean,B,order=10)
}

