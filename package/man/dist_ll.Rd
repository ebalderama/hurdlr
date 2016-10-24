\name{dist_ll}

\alias{dist_ll}

\title{Data Distributional Likelihood for Hurdle Model Count Data Regression}

\description{
\code{dist_ll} is the data likelihood fuction for hurdle model regression using \code{\link{hurdle}}.
}

\usage{
dist_ll(y, hurd = Inf, lam = NULL, size = 1, mu = NULL, xi = NULL, sigma = NULL,
        dist = c("poisson", "nb", "lognormal", "gpd", "none"),
        g.x = F, log = T)
}

\arguments{

\item{y}{
numeric response vector.
}

\item{hurd}{
numeric threshold for 'extreme' observations of two-hurdle models. \code{Inf} for one-hurdle models.
}

\item{lam}{
current value for the poisson likelihood lambda parameter.
}

\item{size}{
size parameter for negative binomial likelihood distributions.
}

\item{mu}{
current value for the negative binomial or log normal likelihood mu parameter.
}

\item{xi}{
current value for the generalized pareto likelihood xi parameter.
}

\item{sigma}{
current value for the generalized pareto likelihood sigma parameter.
}

\item{dist}{
character specification of response distribution.
}

\item{log}{
logical operator. if \code{TRUE}, probabilities p are given as log(p).
}

\item{g.x}{
logical operator. \code{TRUE} if operating within the third component of the likelihood function (the likelihood of 'extreme' observations).
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
