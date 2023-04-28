\name{summary.boodd}
\alias{summary.boodd}

\title{Summary for objects of class boodd}

\description{The bootstrap functions return an object of class \code{boodd}
containing the generated data. The generic function \code{summary} may be
applied to these objects.}

\usage{
\method{summary}{boodd}(object,\dots)
}

\arguments{
\item{object}{an object of class \code{boodd}.}
\item{\dots}{additional arguments.}
}

\details{The function \code{summary.boodd} displays basic information about
the values of the samples computed by the bootsrap functions.
}

\value{The function \code{summary.boodd} returns the \code{kind} attribute
of the function that produced the \code{boodd} object and a (possibly
multicolumn) table with the quartiles, mean, min and max of the computed values of the statistic.
The number of columns is the size of the return value of the function \code{func}.
}

\references{
Efron, B., Tibshirani, R. (1993). \emph{An Introduction to the
Bootstrap}, Chapman and Hall.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.
}


\seealso{\code{\link{confint.boodd}},
\code{\link{plot.boodd}},
\code{\link{boodd-class}}.}


\examples{
B <- 299
x <- round(rnorm(15),3)
boo1 <- boots(x,mean,B)
summary(boo1)

# bootstrap of several statistics 
mv <- function(data) {c(mean(data),var(data))} # compute both mean and variance
boo2 <- boots(x,mv,B)
summary(boo2)

}

