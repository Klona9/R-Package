\name{jackFunc}
\alias{jackFunc}
\alias{jackFuncBlock}

\title{Create a function that calculates the statistic and the jackknife
variance.}

\description{
	Create a vector-valued function that calculates both the statistics
	defined by 'func' and the estimated jackknife variance.
}

\usage{
	jackFunc(func,...)
	jackFuncBlock(func,blen=NULL,...)
}

\arguments{
	\item{func}{the function to apply to each sample.}
	\item{blen}{(integer) block length.}
	\item{...}{optional additional arguments for the \code{func} function.}
}

\details{
	The \code{jackFunc} function creates a vector-valued function that
	calculates both the statistics defined by 'func' and the estimated
	jackknife variance. This function can then be used in generic bootstrap
	procedure to construct bootstrap-t confidence intervals. The
	\code{jackFuncBlock} function is similar but uses the estimated
	jackknife variance based on non-overlapping blocks (calculated with the
	\code{jackVarBlock} function). If the \code{blen} argument is not
	specified, the created function uses, as a default,
	\code{floor(n^(1/3))} where \code{n} is the length of the vector it is
	applied to. }

\value{
\code{jackFunc} and \code{jackFuncBlock} return a function object that takes a single vector argument.
}

\author{
  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{

Efron, B. (1979). Bootstrap methods: an other look at the jackknife,
\emph{Ann. Statist.}, 7, 1-26.

Gray, H., Schucany, W. and Watkins, T. (1972). \emph{The Generalized Jackknife Statistics}. Marcel Dekker, New-York.

Quenouille, M.H. (1949). Approximate tests of correlation in time-series,
\emph{J. Roy. Statist. Soc., Ser. B}, 11, 68-84.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{regenboot}},
\code{\link{blockboot}},
\code{\link{fieldboot}}.}

\examples{
# Create a function to compute both estimation and variance of the empirical skewness
func <- function(x) {mean((x-mean(x))^3)/(mean((x-mean(x))^2)^(3/2))}
x <- arima.sim(list(order=c(1,0,4),ar=0.5,ma=c(0.7,0.4,-0.3,-0.1)),n=101)
jf <- jackFunc(func)
boo1 <- boots(x,jf,299)
# Compute all confidence intervals including bootstrap-t and symmetric
# bootstrap-t methods
confint(boo1,method="all")

jfb <- jackFuncBlock(func,blen=5)
boo2 <- boots(x,jfb,299)
}

