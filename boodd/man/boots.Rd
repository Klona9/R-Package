\name{boots}
\alias{boots}

\title{Bootstrap for the iid case}

\description{Bootstrap for the iid case as described in chapter 2 of [3]. See also [1] and [2] for details.}

\usage{
	boots(x,func,B,smooth=FALSE,moonsize=NULL,mreplace=TRUE,...)
}

\arguments{
\item{x}{a vector or a matrix.}
\item{func}{the function to apply to each sample.}
\item{B}{the number of bootstrap samples.}
\item{smooth}{(logical) specify smooth bootstrap.}
\item{moonsize}{size for 'm out of n' bootstrap.}
\item{mreplace}{(logical) TRUE if bootstrap is done with replacement in 'm out of n'.}
\item{...}{optional additional arguments for the \code{func} function.}
}

\details{The function \code{boots} performs a naive bootstrap in the iid
case. The \code{func} argument must be a function whose first argument is a
vector and which returns either a single value or a vector.

The \code{x} argument can be a vector or a matrix. In the case of a matrix,
the \emph{rows} of the matrix are bootstrapped.

The \code{moonsize} and \code{mreplace} arguments concern \emph{m out of n}
bootstrap (aka \emph{moon bootstrap}). The \code{moonsize} argument is an
integer less than the length of \code{x} (or the number of rows if \code{x}
is a matrix). The \code{mreplace} argument is a logical that indicates
whether the bootstrap samples are drawn with or without replacement. }

\value{An object of class \code{boodd} containing either a vector or a
matrix, depending on whether the \code{func} function returns a single
value or a vector.
If the \code{func} function returns a vector of size \code{n}, then
\code{boots} returns a matrix of size \code{B x n}.
}

\author{
  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{

[1] B. Efron, R. Tibshirani(1993). \emph{An Introduction to the
Bootstrap}, Chapman and Hall.

[2] A. C. Davison, D. Hinkley (1997). \emph{Bootstrap Methods and Their Application}, Cambridge Series in Statistical and Probabilistic Mathematics.

[3] P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{plot.boodd}},
\code{\link{confint.boodd}}.}

\examples{
B <- 299
n <- 200
x <- rnorm(n)
boo1 <- boots(x,mean,B)
summary(boo1)
plot(boo1)
confint(boo1)
confint(boo1,method="bperc")

# bootstrap of several statistics 
mv <- function(data) {c(mean(data),var(data))} # compute both mean and variance
boo2 <- boots(x,mv,B)
# Compute both percentile and t-percentile confidence intervals when variance is bootstrapped
confint(boo2,method="all")

# Naive Bootstrap of all the output parameters of lm (linear regression) function
sigma <- 0.2
y <- x+rnorm(n)
data <- as.matrix(data.frame(x,y))
# the function of interest is here the coeff of the linear model
nlm <- function(dat){lm(dat[,2]~dat[,1])$coefficient}
boo3 <- boots(data,nlm,B)

# Smoothed bootstrap for quantiles 
boo4 <- boots(x,median,B) # without smoothing
plot(boo4)
boo5 <- boots(x,median,B,smooth=TRUE) # with smoothing using a cross-validation estimator of the window
plot(boo5)

# Moon bootstrap 
n <- 10000
x <- rnorm(n)
boo6 <- boots(x,max,B) # i.i.d bootstrap is not consistant for the max
boo7 <- boots(x,max,B,moonsize=sqrt(n)) # a reasonable approximation with moon bootstrap
boo8 <- boots(x,max,B,moonsize=sqrt(n),mreplace=TRUE) # quite similar with or without 
# due to the fact that the probability to observe ties in a bootstrap sample  when moonsize=sqrt(n) is very small. 
}

\keyword{bootstrap}
