\name{boodd}
\alias{boodd}
\alias{boodd-package}

\title{Package overview: bootstrap for dependent data}
\description{
  The package \code{boodd} contains functions, datasets and examples to
  accompany the text book \emph{Bootstrap for Dependent Data} by Patrice Bertail and Anna Dudek.
}

\details{
  \tabular{ll}{
    Version: \tab \Sexpr[stage=render]{library(utils); packageVersion("boodd")}\cr
    Date: \tab \Sexpr[stage=build]{format(Sys.time(),"\%Y-\%m-\%d")}\cr
    License: \tab \Sexpr[stage=build]{library(utils); packageDescription("boodd")$License} \cr
  }

  A list of functions:

\code{\link{boots}} -  Bootstrap in the i.i.d. case. 

\code{\link{bootsemi}} - Semi-parametric bootstrap for time series. 

\code{\link{blockboot}} - Block bootstrap(s) for stationary time series. 

\code{\link{regenboot}} - Regenerative Bootstrap. 

  }

\section{Author}{
  Patrice Bertail \cr
  \email{patrice.bertail@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr 
  Chaire Big Data, TéléComParis-Tech \cr \cr

  Bernard Desgraupes \cr
  \email{bernard.desgraupes@parisnanterre.fr}\cr
  University of Paris Ouest - Nanterre - Lab Modal'X \cr \cr

  Anna Dudek\cr
  Department of Applied Mathematics \cr
  AGH University of Science and Technology - Krakow, Poland \cr \cr
}

\references{P. Bertail, B. Desgraupes and A. Dudek (2023), \emph{Bootstrap for Dependent Data with "R package"}, Springer, N-Y.}

\keyword{package}
