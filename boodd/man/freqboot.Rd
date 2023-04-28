\name{freqboot}
\alias{freqboot}

\title{Frequency Domain Bootstrap}

\description{
	This function implements the Frequency Domain Bootstrap (FDB).
}

\usage{
	freqboot(x,XI,g,B,kernel="normal",bandwidth)
}

\arguments{
	\item{x}{a vector or time series.}
	\item{XI}{a list of functions defined on the interval \eqn{[0,\pi]}.}
	\item{g}{a numeric function taking \code{length(XI)} arguments.}
	\item{B}{the number of bootstrap samples.}
	\item{kernel}{a character string which determines the smoothing kernel.}
	\item{bandwidth}{the kernel bandwidth smoothing parameter.}
}

\details{
The series \code{x} is supposed to be a sample of a real valued zero-mean
stationary time series.

The \code{XI} argument is a list of functions \eqn{\xi_i} (i=1,\dots,p)
used to define a linear functional of the form 
\deqn{A(\xi,f)=(  \int_{0}^{\pi}\xi_{1}(\omega)f(\omega)d\omega,\int_{0}^{\pi
}\xi_{2}(\omega)f(\omega)d\omega,...,\int_{0}^{\pi}\xi_{p}(\omega
)f(\omega)d\omega)
}

The \code{g} argument is a numeric function with \emph{p} arguments so that the statistic
computed by bootstrap is \eqn{T(f)=g(A(\xi,f))}.

An estimate of the spectral density is obtained by smoothing the periodograms of the series. The 
\code{kernel} argument specifies the smoothing kernel. The possible values
are (the names can be abbreviated):
\describe{
	\item{\emph{"normal"}}{the Gaussian density function (the default). }
	\item{\emph{"epanechnikov"}}{the centered beta(2,2) density.}
	\item{\emph{"box"}}{the uniform density on [-1,1].}
}

If the \code{bandwidth} argument is not specified, an optimal value is computed.
Refer to the package for the optimal value.
}

\value{
\code{freqboot} returns an object of class \code{boodd} (see \link{class.boodd}).
}

\author{
  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{
P. Bertail and A.E. Dudek (2021). Consistency of the Frequency Domain
Bootstrap for differentiable functionals, \emph{Electron. J. Statist.}, 15(1), 1-36.

Lahiri, S.N. (2003). \emph{Resampling Methods for Dependent Data}. Springer,
New York.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{aidedboot}}}

\examples{
set.seed(123)
n <- 120
x <- arima.sim(list(order=c(1,0,0),ar=0.7),n=n)
B <- 299
one <- function(x) {1}
XI <- list(cos,one)
g <- function(x,y) {return(x/y)}
# This gives an estimate for the autocorrelation of order 1
freqboot(x,XI,g,B,"normal")
}

\keyword{frequency domain, fft, periodogram}
