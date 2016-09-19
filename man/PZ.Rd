\name{PZ}

\alias{PZ}

\title{Zero Count Probability Likelihood}

\description{
\code{PZ} is used to calculate the likelihood of a zero-value count observation in a single or double-hurdle regression model.
}

\usage{
PZ(p, log = T)
}

\arguments{

\item{p}{
vector of zero-count probabilities.
}
  
\item{log}{
logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
