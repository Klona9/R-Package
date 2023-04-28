\name{jackVar}
\alias{jackVar}
\alias{jackVarBlock}
\alias{jackVarRegen}
\alias{jackVarField}

\title{Jackknife Variance for statistics based on i.i.d data}

\description{
	Estimate the variance of a statistic applied to a vector or a matrix
	using a jackknife procedure.
}

\usage{
	jackVar(x,func,...)
	jackVarBlock(x,func,blen,...)
	jackVarRegen(x,func,...,atom,small=NULL,s=median(x))
	jackVarField(arr,func,blocklens,...)
}

\arguments{
	\item{x}{a vector or a matrix.}
	\item{func}{the function to apply to each sample.}
	\item{arr}{a data array.}
	\item{blocklens}{an integer or a vector of integers for the block lengths.}
	\item{atom}{atom used to cut finite states Markov chains.}
	\item{s}{a real number which is the center for the small set.}
	\item{small}{an object of class \code{smallEnsemble} (see \code{\link{regenboot}}).}
	\item{...}{optional additional arguments for the \code{func} function.}
}

\details{
	If \code{x} is a vector of length \emph{n} or a matrix with \emph{n} rows,
	for all \emph{i} in \code{[[1,n]]}, the statistic \code{func} is calculated on 
	\code{x} with its \emph{i}-th element removed:
	\deqn{T_{n-1}^{i}=func(x[-i])}
	
	Then, for \eqn{T_n=func(x)}:
	\deqn{J^{i}=nT_n - (n-1) T_{n-1}^{i},}
	\deqn{\bar{J}=\frac{1}{n}\sum J^{i},}
	\deqn{V=\frac{1}{n(n-p)}\sum_{i=1}^n (J^{i}-\bar{J}) {}^t(J^{i}-\bar{J}),}
    where \emph{p} is the size of the vector returned by the function \code{func}.

	The \code{jackVarBlock} function computes the same kind of estimator
	but using blocks of observations of length \code{blen}. The sample
	\code{x} is split into
	\deqn{T_{n-1}^{i}=func(x[-B_i]).}

	Then
	\deqn{J^{i}=nT_n - (n-blen) T_{n-1}^{i},}
	\deqn{\bar{J}=\frac{1}{bnum}\sum J^{i},}
	\deqn{V=\frac{1}{n\ bnum}\sum_{i=1}^bnum (J^{i}-\bar{J})^t(J^{i}-\bar{J}).}
	
	The \code{jackVarField} function is similar to \code{jackVarBlock} but
	applied to random fields with mutidimensional blocks.

	The \code{jackVarRegen} function is similar to \code{jackVarBlock} in
	the case of regenerative statistics. It handles variable length blocks
	and supports both the finite states and homogeneous Markov chains. See
	details in the documentation of the \code{\link{regenboot}} and 
	\code{\link{smallEnsemble}} functions.
}	

\value{
\code{jackVar}, \code{jackVarBlock} and \code{jackVarRegen} return a scalar or a variance-covariance matrix depending
on whether the function \code{func} is univariate or multivariate. If the
function \code{func} returns a vector of length \code{p}, the
variance-covariance matrix has size \code{p x p}.
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
set.seed(1)
x <- rnorm(101)
# Function returning a single value
func <- function(x) {sum(abs(x))}
jackVar(x,func)
# Function returning a vector with two elements
mfunc <- function(x) {c(mean(abs(x)),sd(x))}
jackVar(x,mfunc)
# Function with additional argument
funca <- function(x,p) {sum(abs(x)^p)}
jackVar(x,funca,3)

x <- arima.sim(list(order=c(1,0,4),ar=0.5,ma=c(0.7,0.4,-0.3,-0.1)),n=101)
# Jackknife variance estimator of 'func' with blocks of length blen
blen <- 10
V1 <- jackVarBlock(x,func,blen)
V2 <- jackVarBlock(x,mfunc,blen)
V3 <- jackVarBlock(x,funca,blen,2)

# jackVarRegen function in the case of finite state Markov chains
acgt <- c("A","C","G","T")
probs <- c(.3,.1,.3,.3)
n <- 100
atom <- "A"
set.seed(1)
y <- sample(acgt,n,prob=probs,repl=TRUE)
propAtom <- function(x) {
    tbl <- as.vector(table(x))
    prop <- tbl[1]/length(x)
    return(prop)
}
jackVarRegen(y,propAtom,atom=atom)

}

