\name{qVar}
\alias{qVar}

\title{Variance estimator for quantile}

\description{
	This is a kernel based estimator of the variance for quantile.
}

\usage{
	qVar(x,alpha,hn=NULL,kernel=c("gaussian","epanechnikov","rectangular"))
}

\arguments{
	\item{x}{a vector.}
	\item{alpha}{a numeric vector of probabilities.}
	\item{hn}{the smoothing bandwith.}
	\item{kernel}{the smoothing kernel to be used. Possible values are 
		\code{"gaussian"}, \code{"epanechnikov"}, or \code{"rectangular"}
		and may be abbreviated. The default is \code{"gaussian"}.}
}

\details{
Consider the cumulative distribution function \eqn{F} of \eqn{P}, which is
supposed to be continuous and differentiable, with a density \eqn{f}. Denote a
quantile \eqn{T(P)=F^{-1}(\alpha)} for some \eqn{\alpha\in(0,1)}. Assume that
\deqn{
f(F^{-1}(\alpha)) \neq 0.
}
and define the quantile as the unique solution of the equation
\deqn{
F(T(P))=\alpha.
}
Put \eqn{T(P)=F^{-1}(\alpha)} and \eqn{T(P_{n})=F_{n}^{-1}(a)}. In that case the variance estimator is given by
\eqn{S(P)^{2}=\alpha(1-\alpha) / f(F^{-1}(\alpha))^{2}}.

A kernel density estimator is defined by
\deqn{
\hat{f}_{n}^{h_{n,1}}(x) = \frac{1}{nh_{n}}\sum_{i=1}^{n} k(
\frac{x-X_{i}}{h_{n,1}}),
}
where \eqn{h_{n,1}} is the window or smoothing parameter.

Then the estimator of the variance is given by
\deqn{
S_{n,h_{n,1}}^{2}=\frac{\alpha(1-\alpha)}{\hat{f}_{n}^{h_{n,1}}(F_{n}^{-1}(\alpha))^{2}},
}
where the smoothing parameter is such that
\deqn{
nh_{n,1}^{3} / log(n)^{2}  ->  \infty
}
and
\deqn{
nh_{n,1}^{2r+1-\tau}  ->  0
}
for some \eqn{\tau \in (0,1)}.

If the bandwith argument \code{hn} is not specified, a default value of
order \eqn{n^(-1/3)} is computed by unbiased cross-validation. 
}

\value{
\code{qVar} returns a vector of variance estimates of the same length as the argument \code{alpha}.
}

\author{
  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{
Chapter 8 of R.J. Serfling. \emph{Approximation Theorems of Mathematical Statistics}, Wiley Series in Probability and Statistics.

Chapter 3 of P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\examples{
set.seed(123)
x <- rnorm(101)
qVar(x,seq(0, 1, 0.25))
qVar(x,0.25,kernel="epanechnikov")
}

