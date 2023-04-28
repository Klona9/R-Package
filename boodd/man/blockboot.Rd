\name{blockboot}
\alias{blockboot}

\title{Block Bootstrap}

\description{
The function applies block bootstrap methods to a time series.

This function allows the following block bootstrap methods to be
used: the Moving Block Bootstrap (\cite{kunsch} and \cite{LiuSi92}),
the Nonoverlapping Block Bootstrap (\cite{Carlstein1986}), the Circular
Block Bootstrap (\cite{PolRom92}), the Stationary Bootstrap
(\cite{PolRom94}) and the Generalized Seasonal Block Bootstrap
(\cite{DLPP2014}).

}

\usage{
	blockboot(x,func,B,length.block,...,
	     method=c("movingblock","nonoverlapping","circular","stationary","seasonal"),
		 period)
}

\arguments{
	\item{x}{a time series.}
	\item{func}{the function to apply to each sample.}
	\item{B}{the number of bootstrap samples.}
	\item{length.block}{length of blocks.}
	\item{...}{optional additional arguments for the \code{func} function.}
	\item{method}{the block bootstrap method.}
	\item{period}{an integer specifying the length of the period.}
}

\details{
The possible values of the \code{method} argument are: 
\emph{"movingblock"}, \emph{"nonoverlapping"},
\emph{"circular"}, \emph{"stationary"} or \emph{"seasonal"}.
If it is not specified, the default method is \emph{"movingblock"}.
Method names may be abbreviated.

}

\value{
\code{blockboot} returns an object of class \code{boodd}. 
}

\author{
	Bernard Desgraupes \cr
	\email{bernard.desgraupes@parisnanterre.fr}\cr
	University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}

\references{

[kunsch] Künsch, H. (1989). The jackknife and the bootstrap for
general stationary observations. \emph{Ann. Statist.}, 17, 1217-1241.

[LiuSi92] Liu, R. and Singh, K. (1992). Moving block jackknife and
bootstrap capture weak dependence. Exploring the Limits of Bootstrap, Wiley
Ser. Probab. Math. Statist. Probab. Math. Statist. Wiley, New York, pp 225-248.

[PolRom94] Politis, D.N. and Romano, J.P. (1994). The stationary
bootstrap. \emph{J. Amer. Statist. Assoc.}, 89, 1303--1313.

[PolRom92] 
Politis, D.N. and Romano, J.P. (1992).
A circular block-resampling procedure for stationary data.
\emph{Exploring the Limits of Bootstrap},
Wiley Ser. Probab. Math. Statist. Probab. Math. Statist. Wiley, New York, pp 263-270.

[DLPP2014] 
Dudek, A.E., Le\'{s}kow, J., Paparoditis, E. and Politis, D. (2014a). A
generalized block bootstrap for seasonal time series. \emph{J. Time Ser.
Anal.}, 35, 89-114.

[Carlstein1986] Carlstein E. (1986). The use of subseries methods for
estimating the variance of a general statistic from a stationary time series.
\emph{Annals of Statist.}, 14, 1171-1179.


P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{boots}},
\code{\link{bootsemi}},
\code{\link{plot.boodd}},
\code{\link{confint.boodd}},
\code{\link{fieldboot}},
\code{\link{jackVarBlock}}.}

\examples{
B <- 299
data(airquality)
x <- airquality$Wind
n <- length(x)
b <- floor(sqrt(n))
boo1 <- blockboot(x,mean,B,b,method="moving")
plot(boo1,main="MBB")

% The nottem time series
x <- nottem
boo2 <- blockboot(x,mean,B,b,period=12,method="seasonal")
plot(boo2,main="GSBB")
}

\keyword{bootstrap}
