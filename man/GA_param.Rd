\name{GA_param}
\alias{GA_param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GA_param
}
\description{
A list object contains the hyperparameters for GA running.
}
\arguments{
  \item{popsize}{An integer represents the number of individual in each population.}
  \item{Pcrossover}{The probability that the crossover operator applies on two individual chromosomes.}
  \item{Pmutation}{The probability that the mutation operator applies on one individual chromosome.}
  \item{Pchangepoint}{The probability that a changepoint has occurred.}
  \item{minDist}{The minimum length between two adjacent changepoints.}
  \item{mmax}{The maximum possible number of changepoints in the data set.}
  \item{lmax}{The maximum possible length of the chromosome representation.}
  \item{maxgen}{The maximum number of generation that the GA can last.}
  \item{maxconv}{If the overall best fitted value doesn't change after \code{maxconv} consecutive migrations, the island-based GA algorithm stops.}
  \item{option}{A string controls the optimization task. ``cp'' indicates the task is changepoint detection only. ``both'' indicates the task will include both changepoint detection and model order selection.}
  \item{monitoring}{A binary interger 0 or 1, indicating whether print out middle results for each iterations of GA.}
  \item{parallel}{Whether use multiple threads to parallel compute the individual fittness function values.}
  \item{nCore}{An integer represents the number of cores used in parallel computing.}
  \item{tol}{The tolerance level for deciding GA to stop.}
  \item{seed}{An single integer allows function produce reproducible results.}
}
\details{TBD
%%  ~~ If necessary, more details than the description above ~~
}
\references{TBD
%% ~put references to the literature/web site here ~
}
\author{Mo Li}
\examples{
# time series length
Ts = 1000

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
