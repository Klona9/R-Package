\name{class.boodd}
\alias{class.boodd}

\title{Objects of class boodd}

\description{The bootstrap functions return an object of class \code{boodd}
containing the generated data. Some generic functions may be applied to
these objects. See their documentation for more details.}


\details{The functions \code{plot}, \code{confint} and
\code{summary} can be applied directly to all the objects of class
\code{boodd}.}

\value{An object of class \code{boodd} is a list containing at least two components :

\code{$s} is a vector or matrix of the values of the statistic for all the bootstrap samples.
\code{$Tn} is the value of the statistic for the initial series.

}

\references{
Efron, B., Tibshirani, R. (1993). \emph{An Introduction to the
Bootstrap}, Chapman and Hall.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.
}


\seealso{\code{\link{confint.boodd}},
\code{\link{plot.boodd}},
\code{\link{summary.boodd}},
\code{\link{boots}},
\code{\link{blockboot}},
\code{\link{regenboot}},
\code{\link{bootsemi}}.}


\examples{
B <- 299
n <- 200
x <- rnorm(n)
boo1 <- boots(x,mean,B)
print(boo1)
names(boo1)

# bootstrap of several statistics 
mv <- function(data) {c(mean(data),var(data))}
boo2 <- boots(x,mv,B)
print(boo2)

}

