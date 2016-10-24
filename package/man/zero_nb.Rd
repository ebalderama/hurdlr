\name{zero_nb}

\alias{zero_nb}

\title{Zero-Inflated Negative Binomial Regression Model}

\description{
\code{zero_nb} is used to fit zero-inflated negative binomial regression models to count data via Bayesian inference.
}
\usage{
zero_nb(y, x, size, a = 1, b = 1, mu.start = 1,
        beta.prior.mean = 0, beta.prior.sd = 1,
        iters = 1000, burn = 500, nthin = 1,
        plots = T, progress.bar = T)
}
%- maybe also 'usage' for other objects documented here.

\arguments{

\item{y}{
numeric response vector.
}

\item{x}{
numeric predictor matrix.
}

\item{size}{
size parameter for negative binomial likelihood distributions.
}

\item{a}{
shape parameter for gamma prior distributions.
}

\item{b}{
rate parameter for gamma prior distributions.
}

\item{mu.start}{
initial value for mu parameter.
}

\item{beta.prior.mean}{
mu parameter for normal prior distributions.
}

\item{beta.prior.sd}{
standard deviation for normal prior distributions.
}

\item{iters}{
number of iterations for the Markov chain to run.
}

\item{burn}{
numeric burn-in length.
}

\item{nthin}{
numeric thinning rate.
}

\item{plots}{
logical operator. \code{TRUE} to output plots.
}

\item{progress.bar}{
logical operator. \code{TRUE} to print progress bar.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
