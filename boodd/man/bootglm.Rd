\name{bootglm}
\alias{bootglm}

\title{Bootstrap for Generalized Linear Model}

\description{
Parametric bootstrap for generalized linear model.
}

\usage{
	bootglm(model,data,func,B,...)
}

\arguments{
	\item{model}{an object of class \code{lm} or \code{glm}}
	\item{data}{the dataframe used to fit the model}
	\item{B}{the number of bootstrap samples.}
	\item{func}{the function to apply to each sample.}
	\item{...}{optional additional arguments for the \code{func} function}
}

\details{
The parametric bootstrap simply consists in resampling data in the 
model with estimated parameters. \code{bootglm} uses this principle 
for generalized linear models conditionally to the explanatory 
variables (see \cite{beran97})
	
	The \code{model} argument must be a model fitted with the \code{\link{lm}} or \code{\link{glm}} functions.
}

\value{
\code{bootglm} returns an object of class \code{boodd} (see \link{class.boodd}).
}

\author{
	Bernard Desgraupes \cr
	\email{bernard.desgraupes@parisnanterre.fr}\cr
	University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr
}


\references{
[beran97] Beran, R. (1997). \emph{Diagnosing Bootstrap Success}, Annals of 
the Institute of Statistical Mathematics 49, 1-24.

P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\seealso{\code{\link{bootsemi}}.}

\examples{
B <- 299

# Family binomial
x <- runif(100)
e <- 0.5*rnorm(100)
y <- (x + e>0)
data <- data.frame(x,y)
glm_probit <- glm(y ~ x, family=binomial(link="probit"),data=data) 
coeff <- function(data){gg=glm(y ~ x, family=binomial(link="probit"),data=data)$coeff}
  
boo1 <- bootglm(glm_probit,data,coeff,B) # parametric bootstrap of the coeff of a probit model
plot(boo1,main=c("Boostrap of GLM : probit model","Coefficient"))

# coeffv : a function to return coeff and variance of coefficients 
coeffv <- function(data){
    gg <- glm(y~ x, family=binomial(link="probit"),data=data) 
    var <- diag(summary(gg)$cov.u)
    c(gg$coeff,var)
} 
# Parametric bootstrap of all coeff and variances
boo2 <- bootglm(glm_probit,data,coeffv,B)
# Construct all type of confidence intervals including bootstrap-t 
# and symmetric bootstrap-t 
confint(boo2,method="all")

# Family Poisson (Poisson regression)
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
data <- data.frame(treatment,outcome,counts)
gl <- glm(counts ~ outcome + treatment,family=poisson())
meancounts <- function(data) {mean(data$counts)}
boo3 <- bootglm(gl,data,meancounts,B)

}

\keyword{generalized linear model}