\name{get.var}
\alias{get.var}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Get named variable or evaluate expression from list or data.frame}
\description{ This routine takes a text string and a data frame or list. It first sees if the 
string is the name of a variable in the data frame/ list. If it is then the value of this variable is returned. 
Otherwise the routine tries to evaluate the expression within the data.frame/list (but nowhere else) and if 
successful returns the result. If neither step works then \code{NULL} is returned. The routine is useful for
processing gam formulae.}

\usage{ get.var(txt,data)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{txt}{a text string which is either the name of a variable in \code{data} or when 
parsed is an expression that can be evaluated in \code{data}. It can also be neither in which case the
function returns \code{NULL}.}
\item{data}{A data frame or list.} 
}

\value{The evaluated variable, or \code{NULL}}

\references{
\url{http://www.stats.gla.ac.uk/~simon/}
}
\author{ Simon N. Wood \email{simon@stats.gla.ac.uk} } 

\seealso{  \code{\link{gam} } }

\examples{
y <- 1:4;dat<-data.frame(x=5:10)
get.var("x",dat)
get.var("y",dat)
get.var("x==6",dat)
}
\keyword{models} \keyword{smooth} \keyword{regression} %-- one or more ..





