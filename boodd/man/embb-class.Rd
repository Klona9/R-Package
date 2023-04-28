\name{embb}
\alias{embb}
\alias{seasonalMean.embb}
\alias{seasonalVar.embb}
\alias{seasonalACF.embb}
\alias{meanCoeff.embb}
\alias{acfCoeff.embb}

\title{Characteristics for Extended Moving Block Bootstrap class}

\description{
The class of functions that uses the sample obtained with the Extension of
Moving Block Bootstrap (EMBB) method or its circular version (CEMBB). The
functions calculate the seasonal means, seasonal variances and seasonal
autocovariances when a PC time series with period length \emph{d} is considered. For both
PC and APC time series, the functions calculate the Fourier coefficients of
the mean and autocovariance functions.
}

\usage{
	
\method{seasonalMean}{embb}(x,d,\dots)
\method{seasonalVar}{embb}(x,d,\dots)
\method{seasonalACF}{embb}(x,tau,d,\dots)
\method{meanCoeff}{embb}(x,freq,\dots)
\method{acfCoeff}{embb}(x,tau,freq,\dots)

}

\arguments{
	\item{x}{An object of class \code{embb}.}
	\item{d}{A positive integer which is the period length.}
	\item{tau}{Single lag or vector of lags (integers).}
	\item{freq}{Vector of frequencies.}
	\item{\dots}{Additional arguments.}
}

\details{
These methods apply to objects of class \code{embb} typically obtained with the
\code{\link{embb.sample}} functions.
}

\value{
The \code{seasonalMean} and  \code{seasonalVar} functions return a vector of length \code{d}.

\code{seasonalACF} returns either a vector of length \code{d} if a single
lag \code{tau} is specified, or a matrix with \code{length(tau)} rows and \code{d} 
columns if \code{tau} is a vector.

\code{meanCoeff} returns a vector of the same length as \code{freq}.

\code{acfCoeff} returns either a vector of length \code{length(freq)} if a single
lag \code{tau} is specified, or a matrix with \code{length(tau)} rows and \code{length(freq)} 
columns if \code{tau} is a vector.
}

\author{
	Bernard Desgraupes \cr
	\email{bernard.desgraupes@parisnanterre.fr}\cr
	University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{
A.E. Dudek (2015). Circular block bootstrap for coefficients of autocovariance function of almost periodically correlated time series, \emph{Metrika}, 78(3),  313-335.

A.E. Dudek (2018). Block bootstrap for periodic characteristics of periodically correlated
time series. \emph{Journal of Nonparametric Statistics}, 30(1),  87-124.

Anna Dudek, Harry Hurd and Wioletta Wojtowicz (2016). perARMA: Periodic Time Series Analysis. 
https://CRAN.R-project.org/package=perARMA
  
P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{embb.sample}},
\code{\link{seasonalMean.default}},
\code{\link{seasonalVar.default}},
\code{\link{seasonalACF.default}}.
}

\examples{
# Generate a PARMA(2,1) sequence. This requires the perARMA package.
library(perARMA)
T = 12
n = 480
p = 2
a = matrix(0,T,p)
q = 1
b = matrix(0,T,q)

a[1,1] = .8
a[2,1] = .3
phia <- ab2phth(a)
phi0 = phia$phi
phi0 = as.matrix(phi0)

b[1,1] = .7
b[2,1] = .6
thetab <- ab2phth(b)
theta0 = thetab$phi
theta0 = as.matrix(theta0)

del0 = matrix(1,T,1)
PARMA21 <- makeparma(n,phi0,theta0,del0)
x <- PARMA21$y
lb <- 41
# Get an embb object
em <- embb.sample(x,lb,method="movingblock")
d <- 12
seasonalMean(em,d)
seasonalVar(em,d)
h <- c(1,3,6)
seasonalACF(em,h,d)
freq <- 2*pi*(0:(d-1))/d
meanCoeff(em,freq)
acfCoeff(em,h,freq)
}

