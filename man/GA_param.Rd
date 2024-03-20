\name{GA_param}
\alias{GA_param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GA_param
}
\description{
A list object contains the hyperparameters for GA running.
}
\usage{
GA_param = list(
  popsize      = 40,
  Pcrossover   = 0.95,
  Pmutation    = 0.1,
  Pchangepoint = 0.06,
  minDist      = 3,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxgen       = 1000,
  maxconv      = 100,
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
}
\arguments{
  \item{popsize}{An integer represents the number of individual in each population.}
  \item{Pcrossover}{}
  \item{Pmutation}{}
  \item{Pchangepoint}{The probability that a changepoint has occurred.}
  \item{minDist}{The minimum length between two adjacent changepoints.}
  \item{mmax}{The maximum possible number of changepoints in the data set.}
  \item{lmax}{The maximum possible length of the chromosome representation.}
  \item{maxgen}{The maximum number of generation that the GA can last.}
  \item{maxconv}{}
  \item{monitoring}{A binary interger 0 or 1, indicating whether print out middle results for each iterations of GA.}
  \item{parallel}{Whether use multiple threads to parallel compute the individual fittness function values.}
  \item{nCore}{An integer represents the number of cores used in parallel computing.}
  \item{tol}{The tolerance level for deciding GA to stop.}
  \item{seed}{An single integer allows function produce reproducible results.}
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
%%  ~~who you are~~
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
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
