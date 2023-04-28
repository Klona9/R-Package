\name{regenboot}
\alias{regenboot}
\alias{findBestEpsilon}
\alias{smallEnsemble}
\alias{fastNadaraya}

\title{Regenerative Bootstrap}

\description{
Perform regenerative bootstrap on a Markov chain in the atomic case or an
approximate regenerative boostrap in the general Harris case. }

\usage{
	regenboot(x,func,B,...,atom,small=NULL,s=median(x))
	findBestEpsilon(x,s=median(x),plotIt=FALSE)
	fastNadaraya(x,h)
	smallEnsemble(s,eps,delta,trans)

}

\arguments{
	\item{x}{a vector representing a Markov chain.}
	\item{func}{the function to apply to each sample.}
	\item{B}{integer representing the number of bootstrap samples.}
	\item{atom}{either a real number or a string denoting the atom used to cut finite state Markov chains.}
	\item{s}{a real number which represents the center of the small ensemble.}
	\item{small}{an object of class \code{smallEnsemble} (see details below).}
	\item{plotIt}{logical. If TRUE, draw plot of the number of regenerations as a function of the small ensemble radius.}
	\item{eps}{a positive real number, radius for the small ensemble.}
	\item{delta}{a positive real number, lower bound of the transition probability on the small set.}
	\item{trans}{a vector representing the transition density of the Markov chain.}
	\item{h}{a positive real number, bandwith for kernel density estimator of the transition.}
    \item{...}{optional additional arguments for the \code{func} function.}
}

\details{
This function \code{regenboot} implements two different kinds of regenerative bootstrap:
\itemize{
	\item \emph{A regenerative atomic bootstrap} used for finite state Markov chains.
	\item \emph{An approximate regenerative bootstrap} used to bootstrap continuous Markov chains 
	based on a given small set of the form \emph{[s-eps,s+eps]} where \emph{s} is the center and \emph{eps} the radius.
}

One must specify either the \code{atom} argument or the \code{small}
argument. In the first case, \code{atom} is the state used to split the
Markov chain into blocks ending with the atom. In the second case,
\code{small} is an object of class \code{smallEnsemble} representing the
small ensemble (see the \emph{Value} section). Such objects are typically
obtained using the \code{findBestEpsilon} function but may also be constructed
manually using the \code{smallEnsemble} function. By default, \code{s}
is the median of the original vector.

The function \code{fastNadaraya} computes the estimated transition
densities \eqn{p_n(X_i,X_{i+1})} of the Markov chain. It is used in
particular by \code{findBestEpsilon}. It is a Nadaraya kernel type
estimator of the transition density, with a bandwith \emph{h} provided by the
user. The optimal \emph{h} is automatically computed inside the function
\code{findBestEpsilon}.

}

\value{
\code{regenboot} returns an object of class \code{boodd} (see \link{class.boodd}).

\code{findBestEpsilon} and \code{smallEnsemble} return an object of class \code{smallEnsemble}
which is a list with the following components:
\itemize{
	\item \emph{s}: the middle of the small ensemble
	\item \emph{epsilon}: the estimated best radius of the ensemble
	\item \emph{delta}: the estimated lower bound of the transition probability
	\item \emph{trans}: the estimated transition densities \eqn{p_n(X_i,X_{i+1})}
}

\code{fastNadaraya} returns a vector of size \code{length(x)-1}.

}

\author{
  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}


\references{
Bertail, P., Clémençon, S. (2006a).\emph{Regenerative Block
Bootstrap for Markov Chains}. Bernoulli,  12(4):689--712.

P. Bertail and S. Clémençon.
\emph{Regeneration-based statistics for Harris recurrent
  Markov chains}, pages 1--54.
Number 187 in Lecture notes in Statistics. Springer.

Radulovi\'{c}, D. Renewal type bootstrap for Markov chains. Test 13, 147-192.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{boots}},
\code{\link{blockboot}},
\code{\link{plot.boodd}},
\code{\link{confint.boodd}}.}

\examples{
B <- 299
n <- 200

# Atomic Boostrap
acgt <- c("A","C","G","T")
probs <- c(.3,.1,.3,.3)
atom <- "C"
set.seed(1)
x <- sample(acgt,n,prob=probs,repl=TRUE)
propAtom <- function(x) {
    tbl <- as.vector(table(x))
    prop <- tbl[3]/length(x)
    return(prop)
}
boo <- regenboot(x,propAtom,B,atom=atom)
plot(boo)

# Approximate regenerative bootstrap with estimated small set
ar <- arima.sim(list(c(1,0,0),ar=0.6),n=500)
# Find the small ensemble with the largest number of regenerations
sm <- findBestEpsilon(ar,s=0,plotIt=TRUE)
# Approximate regenerative bootstrap of the mean
rboo <- regenboot(ar,mean,small=sm,B=999)
# Plot the corresponding bootstrap distribution
plot(rboo)
# Compute the bootstrap percentile confidence interval
confint(rboo)

}

\keyword{regenerative,bootstrap,Markov chains}
