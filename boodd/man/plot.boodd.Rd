\name{plot.boodd}
\alias{plot.boodd}

\title{Plot objects of class boodd}

\description{The bootstrap functions return an object of class \code{boodd}
containing the generated data. The generic function \code{plot} may be applied to
these objects.}

\usage{
\method{plot}{boodd}(x,with.density=TRUE,which,byrow=FALSE,\dots)
}

\arguments{
\item{x}{an object of class \code{boodd}.}
\item{with.density}{logical. If TRUE, estimated density of the bootstrap distribution is plotted.}
\item{which}{which columns of the data to plot.}
\item{byrow}{logical. If TRUE, display the matrix of histograms by row.}
\item{\dots}{additional arguments for the \code{\link{hist}} function.}
}

\details{The function \code{plot.boodd} plots a histogram (or a matrix of
histograms) of the output data with an estimation of the bootstrap density if
the \code{with.density} argument is set to \code{TRUE}. If the generated
data have more than 1 column, the function displays a matrix of histograms.
The \code{which} argument lets you specify which columns of the
\code{boodd} object should be plotted. The \code{plot.boodd} function
cannot display more than 6 columns. The \code{byrow} argument indicates if
the matrix is organized by rows or by columns. The \code{main} and
\code{xlab} additional parameters may be vectors in order to specify different strings
for each histogram.

}

\value{The function \code{plot.boodd} returns an invisible list containing
the output of the \code{\link{hist}} function. In the case of multiple
histograms, it is a list of lists.

}

\references{
Efron, B., Tibshirani, R. (1993). \emph{An Introduction to the
Bootstrap}, Chapman and Hall.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.
}


\seealso{\code{\link{confint.boodd}},
\code{\link{summary.boodd}},
\code{\link{boodd-class}}.}


\examples{
B <- 299
x <- round(rnorm(15),3)
boo1 <- boots(x,mean,B)
plot(boo1)
confint(boo1)
confint(boo1,method="bperc")

# bootstrap of several statistics 
mv <- function(data) {c(mean(data),var(data))} # compute both mean and variance
boo2 <- boots(x,mv,B)
# Compute both percentile and t-percentile confidence intervals when variance is bootstrapped
confint(boo2,method="all")

}

