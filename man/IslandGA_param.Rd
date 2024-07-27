\name{IslandGA_param}
\alias{IslandGA_param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IslandGA_param
}
\description{
A list object contains the hyperparameters for island-based GA running.
}
\arguments{
  \item{subpopsize}{An integer represents the number of individual in each sub-population (island).}
  \item{Islandsize}{The number of subpopulation (island).}
  \item{Pcrossover}{The probability that the crossover operator applies on two individual chromosomes.}
  \item{Pmutation}{The probability that the mutation operator applies on one individual chromosome.}
  \item{Pchangepoint}{The probability that a changepoint has occurred.}
  \item{minDist}{The minimum length between two adjacent changepoints.}
  \item{mmax}{The maximum possible number of changepoints in the data set.}
  \item{lmax}{The maximum possible length of the chromosome representation.}
  \item{maxMig}{The maximum number of migrations. After \code{maxMig} migrations, the island-based GA algorithm stops.}
  \item{maxgen}{The maximum number of generations that each subpopulation (island) has. The migration will apply after \code{maxgen} generations for each subpopulation(island).}
  \item{maxconv}{If the overall best fitted value doesn't change after  \code{maxconv} consecutive migrations, the island-based GA algorithm stops.}
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

IslandGA_param = list(
  subpopsize   = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.1,
  minDist      = 2,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxMig       = 1000,
  maxgen       = 50,
  maxconv      = 100,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = FALSE, ###
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
}
