\name{Predict.matrix}
\alias{Predict.matrix}
\alias{Predict.matrix.tprs.smooth}
\alias{Predict.matrix.cr.smooth}
\alias{Predict.matrix.ts.smooth}
\alias{Predict.matrix.cs.smooth}
\alias{Predict.matrix.cyclic.smooth}
\alias{Predict.matrix.tensor.smooth}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction methods for smooth terms in a GAM}
\description{ Takes \code{smooth} objects produced by \code{smooth.construct} methods and obtains the matrix mapping 
the parameters associated with such a smooth to the predicted values of the smooth at a set of new covariate values.
}

\usage{
Predict.matrix(object,data)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{object}{ is a smooth object produced by a \code{smooth.construct} method function. The object 
contains all the information required to specify the basis for a term of its class, and this information is
used by the appropriate \code{Predict.matrix} function to produce a prediction matrix for new covariate values.
 Further details are given in \code{\link{smooth.construct}}.}
\item{data}{A data frame containing the values of the (named) covariates at which the smooth term is to be 
evaluated.}
}

\value{ A matrix which will map the parameters associated with the smooth to the vector of values of the smooth 
evaluated at the covariate values given in \code{object}.}


\details{ Smooth terms in a GAM formula are turned into smooth specification objects of 
class \code{xx.smooth.spec} during processing of the formula. Each of these objects is
converted to a smooth object using an appropriate \code{smooth.construct} function. The \code{Predict.matrix} 
functions are used to obtain the matrix that will map the parameters associated with a smooth term to
the predicted values for the term at new covariate values.

Note that new smooth classes can be added by writing a new \code{smooth.construct} method function and a 
corresponding \code{\link{Predict.matrix}} method function: see the example code provided for 
\code{\link{smooth.construct}} for details.}


\references{ 

Wood, S.N. (2000)  Modelling and Smoothing Parameter Estimation
with Multiple  Quadratic Penalties. J.R.Statist.Soc.B 62(2):413-428

Wood, S.N. (2003) Thin plate regression splines. J.R.Statist.Soc.B 65(1):95-114

Wood, S.N. (in press) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass.


\url{http://www.stats.gla.ac.uk/~simon/}
}

\author{Simon N. Wood \email{simon@stats.gla.ac.uk}}

\seealso{ \code{\link{gam}},\code{\link{gamm}}, \code{\link{smooth.construct}} }

\examples{# See smooth.construct examples
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...
