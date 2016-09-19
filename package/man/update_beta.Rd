\name{update_beta}

\alias{update_beta}

\title{MCMC Second-Component Parameter Update Function for Hurdle Model Count Data Regression}

\description{
MCMC algorithm for updating the second-component likelihood parameters in hurdle model regression using \code{\link{hurdle}}.
}

\usage{
update_beta(y, x, hurd, dist, like.part,
            beta.prior.mean, beta.prior.sd,
            beta, XB, beta.acc, beta.tune,
            g.x = F)
}

\arguments{

\item{y}{
numeric response vector of observations within the bounds of the second component of the likelihood function, \eqn{y[0 < y & y < hurd]}.
}

\item{x}{
optional numeric predictor matrix for response observations within the bounds of the second component of the likelihood function, \eqn{y[0 < y & y < hurd]}.
}

\item{hurd}{
numeric threshold for 'extreme' observations of two-hurdle models.
}

\item{dist}{
character specification of response distribution for the second component of the likelihood function.
}

\item{like.part}{
numeric vector of the current second-component likelihood values.
}

\item{beta.prior.mean}{
mu parameter for normal prior distributions.
}

\item{beta.prior.sd}{
standard deviation for normal prior distributions.
}

\item{beta}{
numeric matrix of current second-component parameter values.
}

\item{XB}{
\eqn{x*beta[,1]} product matrix for response observations within the bounds of the second component of the likelihood function, \eqn{y[0 < y & y < hurd]}.
}

\item{beta.acc}{
numeric matrix of current MCMC acceptance rates for second-component parameters.
}

\item{beta.tune}{
numeric matrix of current MCMC tuning values for regression coefficient estimation.
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
\code{\link{dist_ll}}
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
