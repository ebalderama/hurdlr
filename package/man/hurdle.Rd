\name{hurdle}

\alias{hurdle}

\title{Hurdle Model Count Data Regression}

\description{
\code{hurdle} is used to fit single or double-hurdle regression models to count data via Bayesian inference.
}

\usage{
hurdle(y, x = NULL, hurd = Inf,
       dist = c("poisson", "nb", "lognormal", "gpd"),
       dist.2 = c("none", "gpd", "poisson", "lognormal", "nb"),
       control = hurdle_control(...),
       iters = 1000, burn = 500, nthin = 1
       plots = T, progress.bar = T)
}

\arguments{

\item{y}{
numeric response vector.
}

\item{x}{
optional numeric predictor matrix.
}

\item{hurd}{
numeric threshold for 'extreme' observations of two-hurdle models. \code{Inf} for one-hurdle models.
}

\item{dist}{
character specification of response distribution.
}

\item{dist.2}{
character specification of response distribution for 'extreme' observations of two-hurdle models.
}

\item{control}{
list of parameters for controlling the fitting process, specified by \code{\link{hurdle_control}}.
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

\seealso{

\code{\link{hurdle_control}}

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
