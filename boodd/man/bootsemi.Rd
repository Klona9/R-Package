\name{bootsemi}
\alias{bootsemi}

\title{Semiparametric Bootstrap}

\description{
The function performs a semiparametric bootstrap for a general
statistics, using a time-series model.}

\usage{
	bootsemi(x,func,B,...,model=c("ARIMA","GARCH"),params,
	         model.fit=NULL,model.sim=NULL)
}

\arguments{
	\item{x}{a time series.}
	\item{func}{the function to apply to each sample.}
	\item{B}{the number of bootstrap samples.}
    \item{...}{optional additional arguments for the \code{func} function.}
	\item{model}{the chosen model to fit the time series.}
	\item{params}{the model parameters (see below).}
	\item{model.fit}{fitting function for generic model (see below).}
	\item{model.sim}{simulation function for generic model(see below).}
}

\details{
The default basic models currently supported are : \emph{ARIMA} and \emph{GARCH}.

The argument \code{params} specifies the parameters of the chosen model. In
the case of \emph{ARIMA}, it is a vector of the form \emph{c(p,q)} or
\emph{c(p,d,q)}. In the case of \emph{GARCH}, it is a vector of the form
\emph{c(q)} or \emph{c(p,q)} corresponding to an \emph{ARCH(q)} or 
\emph{GARCH(p,q)} model respectively.

Alternatively, one can specify two functions in the \code{model.fit}  and
\code{model.sim} arguments. They are used to implement a generic bootstrap. The 
\code{model.fit} function has the following prototype:
\preformatted{
    model.fit(x,params)
}
It receives the \code{params} argument specified in the \code{bootsemi} function.
It should return an object describing the model (typically a list
containing all the necessary components for the model). The
\code{model.sim} function has the following prototype:
\preformatted{
    model.sim(model,innovations,params)
}
The \code{innovations} argument is a resampled vector of centered residuals.
It should build a new trajectory of the original process using the
data contained in the model object provided by the \code{model.fit} function.

The Examples section below shows how this can be done in the case of a
Threshold Autoregressive (\emph{TAR}) process.

}

\value{
\code{bootsemi} returns an object of class \code{boodd} (see \link{class.boodd}).
}

\author{
  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{boots}},
\code{\link{blockboot}},
\code{\link{bootglm}},
\code{\link{plot.boodd}},
\code{\link{confint.boodd}}.}

\examples{
# An ARIMA(2,1) process
library(stats)
B <- 299
n <- 200
x <- arima.sim(model=list(ar=c(0.8,-0.4),ma=c(0.2)),n=n)
boo1 <- bootsemi(x,mean,B,model="ARIMA",params=c(2,1))
plot(boo1)

\dontrun{
# A GARCH(1,1) process (needs the TSA package)
library(TSA)
x <- garch.sim(alpha=c(1.5,0.4),beta=0.2,n=n)
# bootstrap of several statistics 
mv <- function(data) {c(mean(data),var(data))} # compute both mean and variance
boo2 <- bootsemi(x,mv,B,model="GARCH",params=c(1,1))

# A TAR(1,1,1) process using the generic implementation.
# This example needs the TSA package for the tar and tar.sim functions.
library(TSA)
x <- tar.sim(n=n,Phi1=c(0,0.5),Phi2=c(0,-1.8),thd=-1,d=1,p=1,sigma1=1,sigma2=2)$y
# Define the fitting function. The 'params' argument is of the form c(p1, p2, d).
myfit <- function(x,params) {
    res <- tar(x, p1=params[1], p2=params[2], d=params[3], order.select=FALSE)
    return(res)
}
# Define the simulation function. The 'fit' argument is a TAR model fitted by the tar function.
mysim <- function(fit,innov,params) {
    n <- length(fit$y)
    nx <- tar.sim(fit,ntransient=0,n=n)$y
    return(nx)
}
boo3 <- bootsemi(x,mean,B,params=c(1,1,1),model.fit=myfit,model.sim=mysim)
print(boo3)
}
}

\keyword{semiparametric,bootstrap}
