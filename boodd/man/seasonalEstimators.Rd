\name{seasonalMean}
\alias{seasonalMean}
\alias{seasonalVar}
\alias{seasonalACF}
\alias{seasonalMean.default}
\alias{seasonalVar.default}
\alias{seasonalACF.default}
\alias{seasonalMean.ts}
\alias{seasonalVar.ts}
\alias{seasonalACF.ts}

\title{Characteristics of periodically correlated time series}

\description{
	Calculate estimates of the seasonal means, variances and
	autocovariances of a periodically correlated time series.
}

\usage{
	seasonalMean(x,d,\dots)
	seasonalVar(x,d,\dots)
	seasonalACF(x,tau,d,\dots)

\method{seasonalMean}{default}(x,d,\dots)
\method{seasonalVar}{default}(x,d,\dots)
\method{seasonalACF}{default}(x,tau,d,\dots)

\method{seasonalMean}{ts}(x,d=frequency(x),\dots)
\method{seasonalVar}{ts}(x,d=frequency(x),\dots)
\method{seasonalACF}{ts}(x,tau,d=frequency(x),\dots)
}

\arguments{
	\item{x}{A periodically correlated time series or a vector.}
	\item{d}{a positive integer representing period length.}
	\item{tau}{An integer or a vector of integers representing a single lag or a vector of lags (positive integers).}
	\item{\dots}{Additional arguments.}
}

\details{
The functions \code{seasonalMean} and \code{seasonalVar} calculate
estimates of the seasonal means and variances respectively. The function
\code{seasonalACF} calculates an estimator of the autocovariance for the given
lags. Lags should be positive integers.

}

\value{
The \code{seasonalMean} and  \code{seasonalVar} functions return a vector of length \code{d}.

\code{seasonalACF} returns either a vector of length \code{d} if a single
lag \code{tau} is specified, or a matrix with \code{length(tau)} rows and \code{d} 
columns if \code{tau} is a vector.
}

\author{
	Bernard Desgraupes \cr
	\email{bernard.desgraupes@parisnanterre.fr}\cr
	University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{
P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.

Hurd, H.L., Miamee, A.G. (2007). \emph{Periodically Correlated Random Sequences: Spectral. Theory and Practice}. Wiley.
}

\seealso{\code{\link{blockboot}},
\code{\link{seasonalMean.embb}},
\code{\link{seasonalVar.embb}},
\code{\link{seasonalACF.embb}}.
}

\examples{
# Means
seasonalMean(nottem)  # The period is already in the time-series object nottem
# Variances
seasonalVar(nottem)
# Autocovariances
seasonalACF(nottem,0)
seasonalACF(nottem,1)
seasonalACF(nottem,c(1,3,6))
}

\keyword{seasonal, mean, variance, autocovariance}
