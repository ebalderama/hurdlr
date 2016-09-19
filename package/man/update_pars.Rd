\name{update_pars}

\alias{update_pars}

\title{MCMC Third-Component Parameter Update Function for Hurdle Model Count Data Regression}

\description{
MCMC algorithm for updating the third-component likelihood parameters in hurdle model regression using \code{\link{hurdle}}.
}

\usage{
update_pars(y, hurd, dist, like.part,
            a, b, size, lam, mu, xi, sigma,
            lam.acc, mu.acc, xi.acc, sigma.acc,
            lam.tune, mu.tune, xi.tune, sigma.tune,
            g.x = F)
}

\arguments{

\item{y}{
numeric response vector of observations within the bounds of the third component of the likelihood function, \eqn{y[y ??? hurd]}.
}

\item{hurd}{
numeric threshold for 'extreme' observations of two-hurdle models.
}

\item{dist}{
character specification of response distribution for the third component of the likelihood function.
}

\item{like.part}{
numeric vector of the current third-component likelihood values.
}

\item{a}{
shape parameter for gamma prior distributions.
}

\item{b}{
rate parameter for gamma prior distributions.
}

\item{size}{
size parameter for negative binomial likelihood distributions.
}

\item{lam}{
current value for the poisson likelihood lambda parameter.
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

\item{lam.acc}{
current MCMC acceptance rate for the poisson likelihood lambda parameter.
}

\item{mu.acc}{
current MCMC acceptance rate for the negative binomial or log normal likelihood mu parameter.
}

\item{xi.acc}{
current MCMC acceptance rate for the generalized pareto likelihood xi parameter.
}

\item{sigma.acc}{
current MCMC acceptance rate for the generalized pareto likelihood sigma parameter.
}

\item{lam.tune}{
current MCMC tuning value for the poisson likelihood lambda parameter.
}

\item{mu.tune}{
current MCMC tuning value for the negative binomial or log normal likelihood mu parameter.
}

\item{xi.tune}{
current MCMC tuning value for the generalized pareto likelihood xi parameter.
}

\item{sigma.tune}{
current MCMC tuning value for the generalized pareto likelihood sigma parameter.
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
