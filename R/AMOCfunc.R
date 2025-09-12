#' Random population initialization for AMOC problem
#'
#' Randomly generate the individuals' chromosomes (changepoint configurations)
#' to construct the first generation population for the at most one changepoint
#' (AMOC) problem.
#'
#' @param popSize An integer represents the number of individual in each
#' population for GA (or subpopulation for IslandGA).
#' @param prange Default is \code{NULL} for only changepoint detection. If
#' \code{prange} is specified as a list object, which contains the range of
#' each model order parameters for order selection (integers). The number of
#' order parameters must be equal to the length of \code{prange}.
#' @param N The length of time series.
#' @param minDist The minimum length between two adjacent changepoints.
#' @param pchangepoint The probability that a changepoint can occur.
#' @param mmax The maximum possible number of changepoints in the data set.
#' @param lmax The maximum possible length of the chromosome representation.
#' @details
#' Given the possible candidate changepoint location set, each chromosome in the
#' first generation population can be obtained by randomly sampling one location
#' from the candidate set. The first element of every chromosome represent the
#' number of changepoints and the last non-zero element always equal to the
#' length of time series plus one (N+1).
#' @return A matrix that contains each individual's chromosome.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
AMOCpopulation <- function(popSize, prange, N, minDist, pchangepoint, mmax, lmax) {
  tauclc <- floor(0.05 * N):ceiling(0.95 * N)

  pop <- matrix(0, nrow = lmax, ncol = popSize)
  pop[1, ] <- rep(1, popSize)
  pop[2, ] <- sample(tauclc, size = popSize)
  pop[3, ] <- N + 1

  return(pop)
}

#' The parents selection genetic algorithm operator for AMOC problem
#'
#' The genetic algorithm require to select a pair of chromosomes, representing
#' \code{dad} and \code{mom}, for the \code{crossover} operator to
#' produce offspring (individual for next generation). Here, the same linear
#' ranking method in \code{\link{selection_linearrank}} is used to select a pair
#' of chromosomes for \code{dad} and \code{mom} in the at most one changepoint
#' (AMOC) problem. By default, the dad has better fit/smaller fitness function
#' value/larger rank than \code{mom}.
#' @param pop A matrix contains the chromosomes for all individuals. The number of
#' rows is equal to \code{lmax} and the number of columns is equal to the
#' \code{popSize}.
#' @param popFit A vector contains the objective function value (population fit)
#' being associated to each individual chromosome from above.
#' @return A list contains the chromosomes for \code{dad} and \code{mom}.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
AMOCselection <- function(pop, popFit) {
  return(selection_linearrank(pop, popFit))
}

#' Average crossover operator to produce offspring for AMOC problem
#'
#' In this crossover operator designed for AMOC problem, the new child is
#' produced by taking the average of the changepoint locations from \code{dad}
#' and \code{mom} and round to an integer. Note, every chromosome has at most
#' one candidate changepoint location.
#' @param mom Among two selected individuals, \code{mom} represents the selected
#' chromosome representation with lower fitness function value.
#' @param dad Among two selected individuals, \code{dad} represents the selected
#' chromosome representation with larger fitness function value.
#' @param prange The default value is \code{NULL}. If there is no requirement
#' on model order selection, such an auxiliary argument is needed for \code{GA}
#' and \code{IslandGA} functions.
#' @param minDist The minimum length between two adjacent changepoints.
#' @param lmax The maximum possible length of the chromosome representation.
#' @param N The length of time series.
#' @return The child chromosome that produced from \code{mom} and \code{dad} for
#' next generation.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
AMOCcrossover <- function(mom, dad, prange = NULL, minDist, lmax, N) {
  child <- matrix(0, nrow = lmax, 1)
  child[1] <- 1
  child[2] <- round((dad[2] + mom[2]) / 2)
  child[3] <- N + 1

  return(child)
}

#' Jump mutation operator to produce offspring for AMOC problem
#'
#' In a certain probability, the \code{mutation} genetic operator can be applied
#' to generate a new \code{child}. In this AMOC mutation operator, the new child
#' changepoint location can be down via a "jump" method. The child changepoint
#' location will jump \code{minDist} time units either to the left or right to
#' produce the changepoint location for the mutated child. The jump direction is
#' randomly decided with 0.5 probability.
#'
#' @param child The child chromosome resulting from the \code{crossover} genetic
#' operator.
#' @param prange The default value is \code{NULL}. If there is no requirement
#' on model order selection, such an auxiliary argument is needed for \code{GA}
#' and \code{IslandGA} functions.
#' @param minDist The minimum length between two adjacent changepoints in
#' \code{\link{AMOCselection}} operator, which is also the jump magnitude in the
#' \code{AMOCmutation} operator.
#' @param pchangepoint An auxiliary argument is needed for \code{GA}
#' and \code{IslandGA} functions.
#' @param lmax An auxiliary argument is needed for \code{GA} and \code{IslandGA}
#' functions.
#' @param mmax An auxiliary argument is needed for \code{GA} and \code{IslandGA}
#' functions.
#' @param N An auxiliary argument is needed for \code{GA} and \code{IslandGA}
#' functions.
#' @return The resulting child chromosome representation.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
AMOCmutation <- function(child, prange = NULL, minDist, pchangepoint = NULL, lmax = NULL, mmax = NULL, N = NULL) {
  tmptau <- 1

  while (tmptau < floor(0.05 * N) | tmptau > ceiling(0.95 * N)) {
    tmpsign <- sample(x = c(-1, 1), size = 1)
    tmptau <- child[2] + tmpsign * minDist
  }
  child[2] <- tmptau

  return(child)
}
