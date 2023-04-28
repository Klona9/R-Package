\name{embb.sample}
\alias{embb.sample}

\title{EMBB method}

\description{
The function constructs the bootstrap sample using the Extension of
Moving Block Bootstrap (EMBB) method, which is valid for periodically
correlated and almost periodically correlated time series.
}

\usage{
	embb.sample(x,length.block,method=c("movingblock","circular"))
}

\arguments{
	\item{x}{almost periodically correlated time series.}
	\item{length.block}{length of blocks.}
	\item{method}{bootstrap method.}
}

\details{
The argument \code{method} can be set to \emph{"movingblock"} (in the case
of the \emph{EMBB}) or to \emph{"circular"} (in the case of the circular
version of the \emph{EMBB}).
Method names may be abbreviated. 
}

\value{
\code{embb.sample} returns a matrix whose first column is the bootstrapped
sample and second column contains the original time indices of the 
chosen observations.
}

\author{
	Bernard Desgraupes \cr
	\email{bernard.desgraupes@parisnanterre.fr}\cr
	University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{
A.E. Dudek (2015). Circular block bootstrap for coefficients of autocovariance function of almost periodically correlated time series, \emph{Metrika}, 78(3),  313-335.

A.E. Dudek (2018).   Block bootstrap for periodic characteristics of periodically correlated
time series. \emph{Journal of Nonparametric Statistics}, 30(1),  87-124.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.
}

\seealso{\code{\link{embb}},
\code{\link{blockboot}}.}

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
x <- makeparma(n,phi0,theta0,del0)$y
lb <- 41

# Moving blocks method
embb.sample(x,length.block=lb,method="movingblock")

# Circular method
embb.sample(x,length.block=lb,method="circular")
}

