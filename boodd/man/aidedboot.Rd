\name{aidedboot}
\alias{aidedboot}

\title{Aided Frequency Bootstrap}

\description{
The Aided Frequency Bootstrap (AFB) estimates functionals of the spectral
density by bootstrapping the periodograms using the quotient of two
spectral densities.

The idea underlying the Aided Frequency Bootstrap is
importance sampling. It was introduced by Kreiss and Paparoditis (2003) and
allows to better mimic the asymptotic covariance structure of the
periodogram in the bootstrap world. Kreiss and Paparoditis (\cite{KP2003})
considered a spectral density which is easy to estimate (typically based on
a sieve AR representation of the time series), say \eqn{f_{AR}(\omega)}.
Then putting
\eqn{q(\omega)=\frac{f(\omega)}{f_{AR}(\omega)}}, the problem is to
estimate this quantity to generate bootstrap values of the periodogram.

}

\usage{
	aidedboot(x,XI,g,B,order=NULL,kernel="normal",bandwidth)
}

\arguments{
	\item{x}{a vector or time series.}
	\item{XI}{a list of functions defined on the interval \eqn{[0,\pi]}.}
	\item{g}{a numeric function taking \code{length(XI)} arguments.}
	\item{B}{the number of bootstrap samples.}
	\item{order}{(integer) the order of the autoregressive sieve process.}
	\item{kernel}{a character string which determines the smoothing kernel.}
	\item{bandwidth}{the kernel bandwidth smoothing parameter.}
}

\details{
The argument \code{x} is supposed to be a sample of a real valued zero-mean
stationary time series.

The autoregressive sieve process of order \eqn{l=l(n)} is modelled as
\deqn{X_{t}=\sum_{k=1}^{l}\psi_{k}X_{t-k}+\epsilon_{t}
}
with \eqn{E(\epsilon_{t})=0, Var(\epsilon_{t})=\sigma^{2}(l)}.

We estimate functionals of the spectral density \eqn{T(f)} of the form
\deqn{
T(f)=g(A(\xi,f))
}
where \eqn{g} is a third order differentiable function,\
\deqn{A(\xi,f)=(  \int_{0}^{\pi}\xi_{1}(\omega)f(\omega)d\omega,\int_{0}^{\pi
}\xi_{2}(\omega)f(\omega)d\omega,...,\int_{0}^{\pi}\xi_{p}(\omega
)f(\omega)d\omega)
}
and
\deqn{
\xi=(\xi_{1},\dots,\xi_{p}): [0,\pi] -> R^p.
}


If the \code{order} argument is not specified, its default value is 
\eqn{floor(4*(n/log(n))^(1/4))}.

The \code{kernel} argument has the same meaning as with the \code{\link{freqboot}}
function.
}

\value{
\code{aidedboot} returns an object of class \code{boodd} (see \link{class.boodd}).
}

\author{
  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{
[KP2003] Kreiss, J.-P. and Paparoditis, E. (2003). Autoregressive aided
periodogram bootstrap for time series. \emph{Ann. Stat.} \bold{31} 1923--1955.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.
}

\seealso{\code{\link{freqboot}}}


\examples{
n <- 200
x <- arima.sim(list(order=c(4,0,0),ar=c(0.7,0.4,-0.3,-0.1)),n=n)
B <- 299
one <- function(x) {1}
XI <- list(cos,one)
g <- function(x,y) {return(x/y)}
ord <- floor(4*n^0.25/log(n)^0.5)
boo <- aidedboot(x,XI,g,B,order=ord)
}

\keyword{bootstrap, frequency domain, periodogram}
