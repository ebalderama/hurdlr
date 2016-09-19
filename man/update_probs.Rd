\name{update_probs}

\alias{update_probs}

\title{MCMC Probability Update Function for Hurdle Model Count Data Regression}

\description{
MCMC algorithm for updating the likelihood probabilities in hurdle model regression using \code{\link{hurdle}}.
}

\usage{
update_probs(y, x, hurd, p, q,
             beta.prior.mean, beta.prior.sd,
             pZ, pT, pE,
             beta, XB2, XB3, beta.acc, beta.tune)
}

\arguments{

\item{y}{
numeric response vector.
}

\item{x}{
optional numeric predictor matrix.
}

\item{hurd}{
numeric threshold for 'extreme' observations of two-hurdle models.
}

\item{p}{
numeric vector of current 'p' probability parameter values for zero-value observations. 
}

\item{q}{
numeric vector of current 'q' probability parameter values for 'extreme' observations. 
}

\item{beta.prior.mean}{
mu parameter for normal prior distributions.
}

\item{beta.prior.sd}{
standard deviation for normal prior distributions.
}

\item{pZ}{
numeric vector of current 'zero probability' likelihood values.
}

\item{pT}{
numeric vector of current 'typical probability' likelihood values.
}

\item{pE}{
numeric vector of current 'extreme probability' likelihood values.
}

\item{beta}{
numeric matrix of current regression coefficient parameter values.
}

\item{XB2}{
\eqn{x*beta[,2]} product matrix.
}

\item{XB3}{
\eqn{x*beta[,3]} product matrix.
}

\item{beta.acc}{
numeric matrix of current MCMC acceptance rates for regression coefficient parameters.
}

\item{beta.tune}{
numeric matrix of current MCMC tuning values for regression coefficient estimation.
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
