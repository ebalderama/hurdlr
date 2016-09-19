\name{hurdle_control}

\alias{hurdle_control}

\title{Control Parameters for Hurdle Model Count Data Regression}

\description{
Various parameters for fitting control of hurdle model regression using \code{\link{hurdle}}.
}

\usage{
hurdle_control(a = 1, b = 1, size = 1,
               beta.prior.mean = 0, beta.prior.sd = 1000,
               beta.tune = 1, pars.tune = 0.2,
               lam.start = 1, mu.start = 1,
               sigma.start = 1, xi.start = 1)
}

\arguments{

\item{a}{
shape parameter for gamma prior distributions.
}

\item{b}{
rate parameter for gamma prior distributions.
}

\item{size}{
size parameter for negative binomial likelihood distributions.
}

\item{beta.prior.mean}{
mu parameter for normal prior distributions.
}

\item{beta.prior.sd}{
standard deviation for normal prior distributions.
}

\item{beta.tune}{
Markov chain tuning for regression coefficient estimation.
}

\item{pars.tune}{
Markov chain tuning for parameter estimation of 'extreme' observations distribution.
}

\item{lam.start}{
initial value for the poisson likelihood lambda parameter.
}

\item{mu.start}{
initial value for the negative binomial or log normal likelihood mu parameter.
}

\item{sigma.start}{
initial value for the generalized pareto likelihood sigma parameter.
}

\item{xi.start}{
initial value for the generalized pareto likelihood xi parameter.
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
