\name{meanCoeff}
\alias{meanCoeff}
\alias{acfCoeff}
\alias{meanCoeff.default}
\alias{acfCoeff.default}
\alias{meanCoeff.ts}
\alias{acfCoeff.ts}

\title{Estimation of the Fourier coefficients.}

\description{
For both periodically (PC) and almost periodically correlated (APC) data,
the functions calculate the Fourier coefficients of the mean and
autocovariance functions. The function can also be used for bootstrap
samples obtained with the EMBB, CEMBB, GSBB, CGSBB.

}

\usage{
	meanCoeff(x,d,freq=NULL,\dots)
	acfCoeff(x,tau,d,freq=NULL,\dots)

\method{meanCoeff}{default}(x,d,freq=NULL,\dots)
\method{acfCoeff}{default}(x,tau,d,freq=NULL,\dots)

\method{meanCoeff}{ts}(x,d=frequency(x),freq=NULL,\dots)
\method{acfCoeff}{ts}(x,tau,d=frequency(x),freq=NULL,\dots)
}

\arguments{
	\item{x}{Periodically or almost periodically correlated time series.}
	\item{d}{(Integer) Period length. By default, \code{frequency(x)}.}
	\item{tau}{Single lag or vector of lags (integers).}
	\item{freq}{Vector of frequencies.}
	\item{\dots}{Additional arguments.}
}

\details{
If the \code{freq} argument is not specified, the Fourier frequencies are used:
\eqn{2*k*pi/d} for \emph{k=0,1,\dots,d} where \emph{d} is the frequency of the time series.

The \code{meanCoeff} function implements the estimator of the Fourier coefficient of the mean
at frequency \eqn{\gamma}:
\deqn{\widehat{b}(\gamma) = \frac{1}{n}\sum_{t=1}^n X_t e^{-i\gamma t}}{
widehat(b)(gamma) = (1/n)sum(X[t] exp(-i \gamma t), t=1, n)}

The \code{acfCoeff} function implements the estimator of the Fourier
coefficient of the autocovariance for given lag \code{tau} at frequency
\eqn{\lambda}:
\deqn{\widehat{a}(\lambda,\tau) = \frac{1}{n}\sum_{t=1-\min\{\tau,0\}}^{n-\max\{\tau,0\}} (X_{t+\tau}-\widehat{\mu}_n(t+\tau))
 (X_{t}-\widehat{\mu}_n(t))e^{-i\lambda t}
}
}

\value{
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

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{Methods for class \code{embb}: \code{\link{meanCoeff.embb}},
\code{\link{acfCoeff.embb}}.
}

\examples{
# Fourier frequencies for the data nottem (temperatures at Nottingham Castle)
meanCoeff(nottem)
acfCoeff(nottem,5)
# Given frequencies
meanCoeff(nottem,freq=1:6)
acfCoeff(nottem,tau=c(1,2,5),freq=1:6)
}

\keyword{frequency domain, estimators}
