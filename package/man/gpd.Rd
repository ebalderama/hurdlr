\name{Gen Pareto}

\alias{dgpd}
\alias{mgpd}
\alias{pgpd}
\alias{qgpd}
\alias{rgpd}

\title{The Generalized Pareto Distribution}

\description{
Density, distribution function, quantile function and random generation for the Generalized Pareto distribution with parameters \code{mu}, \code{sigma}, and \code{xi}.
}

\usage{
dgpd(x, mu = 0, sigma = 1, xi = 1, log = F)
mgpd(x, mu = 0, sigma = 1, xi = 1, log = F)
pgpd(q, mu = 0, sigma = 1, xi = 1, lower.tail = T)
qgpd(p, mu = 0, sigma = 1, xi = 1, lower.tail = T)
rgpd(n, mu = 0, sigma = 1, xi = 1)
}

\arguments{

\item{x, q}{
vector of quantiles.
}
  
\item{p}{
vector of probabilities.
}
  
\item{n}{
number of random values to return.
}
  
\item{mu}{
location parameter.
}
  
\item{sigma}{
(non-negative) scale parameter.
}
  
\item{xi}{
shape parameter
}
  
\item{log}{
logical; if \code{TRUE}, probabilities p are given as log(p).
}
  
\item{lower.tail}{
logical; if \code{TRUE}, probabilities are \eqn{P[X ??? x]}, otherwise, \eqn{P[X > x]}.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
\code{dgpd} gives the continuous density, \code{mgpd} gives the discrete density, \code{pgpd} gives the distribution function, \code{qgpd} gives the quantile function, and \code{rgpd} generates random deviates.
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
Taylor Trippe <\email{ttrippe@luc.edu}>

Dr. Earvin Balderama <\email{ebalderama@luc.edu}>
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hurdle}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
