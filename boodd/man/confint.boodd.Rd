\name{confint.boodd}
\alias{confint.boodd}

\title{Confidence intervals for objects of class boodd}

\description{The bootstrap functions return an object of class \code{boodd}
containing the generated data. The generic function \code{confint} may be
applied to these objects.}

\usage{
\method{confint}{boodd}(object,alpha=0.05,method=c("perc","bperc","aboot","tboot","tsymboot","all"),recenter,\dots)
}

\arguments{
\item{object}{an object of class \code{boodd}.}
\item{alpha}{error level for confidence interval.}
\item{method}{method used to build the confidence interval. The default is \emph{perc}.}
\item{recenter}{logical. If TRUE, indicates whether to center the intervals on the 
mean value of the bootstrap samples (only for \emph{tboot} or \emph{tsymboot} methods).}
\item{\dots}{additional arguments.}
}

\details{
The function \code{confint.boodd} provides a confidence interval. Several
methods are available (see \cite{EfronTibshirani1993}):
\itemize{
	\item \emph{perc}: percentile
	\item \emph{bperc}: basic percentile
	\item \emph{aboot}: asymptotic bootstrap
	\item \emph{tboot}: bootstrap-t
	\item \emph{tsymboot}: symmetric bootstrap-t
	\item \emph{all}: all the previous methods
}

The \emph{tboot} and \emph{tsymboot} methods require the function 
\code{func} to which the bootstrap method is applied to return an 
even number of values corresponding to estimates of the parameters 
and estimates of their variances.

The \emph{recenter} argument is used only by the \emph{tboot} and
\emph{tsymboot} methods and is \emph{TRUE} by default (unless the
object was obtained by \emph{Circular Block Bootstrap}, in which case
it is \emph{FALSE} by default).
}

\value{
If the \code{method} argument is not "all", the function
\code{confint.boodd} returns a two-column matrix representing the lower and
upper bounds of the interval. Each row of the matrix corresponds to the
variable to which the interval applies. The default value of the
\code{method} argument is "perc". If the \code{method} argument is "all",
the function \code{confint.boodd} returns a list with the confidence
intervals for all supported methods.

}

\references{
[EfronTibshirani1993] Efron, B., Tibshirani, R. (1993). \emph{An Introduction to the
Bootstrap}, Chapman and Hall.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.
}


\seealso{\code{\link{plot.boodd}},
\code{\link{summary.boodd}},
\code{\link{boodd-class}}.}


\examples{
B <- 299
x <- round(rnorm(15),3)
boo1 <- boots(x,mean,B)
confint(boo1)
confint(boo1,method="bperc")

# bootstrap of several statistics 
mv <- function(data) {c(mean(data),var(data))} # compute both mean and variance
boo2 <- boots(x,mv,B)
# Compute both percentile and t-percentile confidence intervals when variance is bootstrapped
confint(boo2,method="all")

}

