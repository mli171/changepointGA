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
  \item{popsize}{An integer represents the number of individuals in each population.}
  \item{Pcrossover}{The probability that the crossover operator applies on two individual chromosomes.}
  \item{Pmutation}{The probability that the mutation operator applies on one individual chromosome.}
  \item{Pchangepoint}{The probability that a changepoint has occurred. User
  could change this probability based on domain knowledge and the time series
  length.}
  \item{minDist}{The minimum length between two adjacent changepoints.}
  \item{mmax}{The maximum possible number of changepoints in the data set.
  For a time series of length 1000, the default value is 499. The suggested
  value should be based on the length of the time series. For instance, if a
  time series has length N, the recommended \code{mmax} should be N/2-1.}
  \item{lmax}{The maximum possible length of the chromosome representation. For
  a time series of length 1000, the default value is 501. The suggested value
  should be based on the length of the time series. For instance, if a time
  series has length N, the recommended \code{lmax} should be 2+N/2-1.}
  \item{maxgen}{The maximum number of generation that the GA can last.}
  \item{maxconv}{If the overall best fitted value doesn't change after \code{maxconv} consecutive migrations, the GA algorithm stops.}
  \item{option}{A string controls the optimization task. ``cp'' indicates the task is changepoint detection only. ``both'' indicates the task will include both changepoint detection and model order selection.}
  \item{monitoring}{A binary interger 0 or 1, indicating whether print out middle results for each iterations of GA.}
  \item{parallel}{Whether use multiple threads to parallel compute the individual fittness function values.}
  \item{nCore}{An integer represents the number of cores used in parallel computing.}
  \item{tol}{The tolerance level for deciding GA to stop.}
  \item{seed}{An single integer allows function produce reproducible results.}
}
\author{Mo Li}
\examples{
# time series length
N = 1000

GA_param = list(
  popsize      = 40,
  Pcrossover   = 0.95,
  Pmutation    = 0.1,
  Pchangepoint = 0.06,
  minDist      = 3,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxgen       = 50000,
  maxconv      = 5000,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
}
