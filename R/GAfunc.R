#-------------------------- Check Convergence
check_conv <- function(a, maxconv, tol) {
  a <- tail(a, maxconv)
  
  if (length(a) < 2L) {
    return(0L)
  }
  if (anyNA(a) || any(!is.finite(a))) {
    return(0L)
  }
  
  as.integer(all(abs(diff(a)) < tol))
}

#' The default mutation operator in genetic algorithm
#'
#' In a certain probability, the \code{mutation} genetic operator can be applied
#' to generate a new \code{child}. By default, the new child selection can be
#' down by the similar individual selection method in population initialization,
#' \code{\link{select_tau}}.
#'
#' @param child The child chromosome resulting from the \code{crossover} genetic
#' operator.
#' @param prange Default is \code{NULL} for only changepoint detection. If
#' \code{prange} is specified as a list object, which contains the range of
#' each model order parameters for order selection (integers). The number of
#' order parameters must be equal to the length of \code{prange}.
#' @param minDist The required minimum distance between two adjacent changepoints.
#' @param pchangepoint The probability of changepoints for every time series.
#' @param lmax The user specified maximum number of changepoints, by default,
#' as \code{N/2 - 1}.
#' @param mmax The user specified maximum length of individual chromosome,
#' by default, as \code{2+N/2 + 1}.
#' @param N The sample size of the time series.
#' @details
#' A function can apply mutation to the produced child with the specified
#' probability \code{pmutation} in \code{cptga} and
#' \code{cptgaisl}. If order selection is not requested
#' (\code{option = "cp"} in \code{cptga} and \code{cptgaisl}), the default
#' \code{\link{mutation}} operator function uses \code{select_tau} to select
#' a completely new individual with a new chromosome as the mutated child.
#' For details, see \code{\link{select_tau}}. If order selection is needed
#' (\code{option = "both"} in \code{cptga} and \code{cptgaisl}), we first
#' decide whether to keep the produced child's model order with a probability
#' of 0.5. If the child's model order is retained, the \code{select_tau}
#' function is used to select a completely new individual with a new chromosome
#' as the mutated child. If a new model order is selected from the candidate
#' model order set, there is a 0.5 probability to either select a completely new
#' individual with new changepoint locations or retain the original child's
#' changepoint locations for the mutated child. Note that the current model
#' orders in the child's chromosome are excluded from the set to avoid redundant
#' objective function evaluation. Finally, the function returns a vector
#' containing the modified chromosomes for mutated \code{child}.
#' @return The resulting child chromosome representation.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
mutation <- function(child, prange = NULL, minDist, pchangepoint, lmax, mmax, N) {
  plen <- length(prange)

  if (plen > 0) {
    childMut <- matrix(0, nrow = lmax, 1)
    a1 <- runif(1)
    if (a1 > 0.5) {
      # 1. order from child
      childMut[2:(plen + 1), ] <- child[2:(plen + 1)]
      # 1.1 cpt from new
      tmpchildMut <- select_tau(
        N = N, prange = NULL, minDist = minDist, pchangepoint = pchangepoint,
        mmax = mmax, lmax = lmax
      )
      childMut[1, ] <- tmpchildMut[1]
      childMut[(plen + 2):(plen + tmpchildMut[1] + 2), ] <- tmpchildMut[2:(tmpchildMut[1] + 2)]
    } else {
      # 2. order from new
      new.prange <- rep(NA, plen)
      for (ii in 1:plen) {
        tmp.prange <- setdiff(prange[[ii]][1]:prange[[ii]][2], child[2 + ii - 1])
        new.prange[ii] <- sample(tmp.prange, 1)
      }
      childMut[2:(plen + 1), ] <- new.prange
      a2 <- runif(1)
      if (a2 > 0.5) {
        # 2.1 cpt from new
        tmpchildMut <- select_tau(
          N = N, prange = NULL, minDist = minDist, pchangepoint = pchangepoint,
          mmax = mmax, lmax = lmax
        )
        childMut[1, ] <- tmpchildMut[1]
        childMut[(plen + 2):(plen + tmpchildMut[1] + 2), ] <- tmpchildMut[2:(tmpchildMut[1] + 2)]
      } else {
        # 2.2 cpt from child
        childMut[1, ] <- child[1]
        childMut[(plen + 2):(plen + child[1] + 2), ] <- child[(plen + 2):(plen + child[1] + 2)]
      }
    }
  } else {
    tmpchildMut <- select_tau(
      N = N, prange = NULL, minDist = minDist, pchangepoint = pchangepoint,
      mmax = mmax, lmax = lmax
    )
    childMut <- tmpchildMut
  }

  return(childMut)
}

new_population_Island <- function(ObjFunc, prange, selection, crossover, mutation, pop, fit, minDist, lmax, mmax, pcrossover, pmutation, pchangepoint, maxgen, N, ...) {
  # This function is used to form new population
  # some inputs ++++++++++++++++++
  #   pop= population
  #   fit= fitness evaluated for population
  #   minDist= minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   mmax= max number of changepoints
  #   pcrossover= prob of crossover
  #   pmutation= prob of mutation
  #   pchangepoint= prob of changepoints for every time series
  #   maxgen= for each subpopulation, after maxgen then apply migration
  #   N= sample size
  #   X_hour= categorical time series
  #   XMat= Design matrix including covariate other than changepoint
  #   penalty= selection criterion to choose
  # outputs ++++++++++++++++++
  #   pop= updated population
  #   fit= updated population fitness
  #   bestfit  = currnt minimum of fitness values
  #   bestchrom = chromosome representation of the individual associated with bestfit

  plen <- length(prange)

  count <- 1
  repeat{
    # indicator for c("crossover", "mutation")
    #     flag[1]=1 indicating no cross-over
    #     flag[2]=1 indicating no mutation
    flag <- rep(0, 2)
    ##### step 2: parents selection
    parents <- selection(pop, fit)
    dad <- parents$dad
    mom <- parents$mom

    ##### step 3: crossover
    a1 <- runif(1)
    if (a1 <= pcrossover) {
      child <- crossover(mom, dad, prange, minDist, lmax, N)
    } else {
      child <- dad
      flag[1] <- 1
    }

    ## step 4-2: mutation
    a2 <- runif(1)
    if (a2 <= pmutation) {
      child <- mutation(child, prange, minDist, pchangepoint, lmax, mmax, N)
    } else {
      flag[2] <- 1
    }

    ## step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum <- flag[1] + flag[2]
    if (flagsum < 2) {
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      fitChild <- do.call(ObjFunc, c(list(child[1:(child[1] + plen + 2)], plen, ...)))
      # fitChild = do.call(ObjFunc, c(list(child[1:(child[1]+plen+2)], plen, XMat, Xt)))
      leastfit <- max(fit) # with largest fitness value

      if (fitChild < leastfit) {
        # indicating child is better than the worst one and replace
        pp <- which.max(fit)
        pop[, pp] <- child
        fit[pp] <- fitChild
      }
      count <- count + 1
    }

    # check: after every maxgen generations, apply migration in GA.Main()
    if (count >= maxgen) {
      break
    }
  }

  return(rbind(fit, pop))
}


#' Birth-death-relocate mutation operator
#'
#' Perform a structured birth-death-relocate (BDR) mutation for changepoint
#' chromosomes. Conditional on mutation, the operator randomly selects one
#' feasible operation from birth, death, and relocation, with equal probability
#' among the feasible operation types.
#'
#' The operator can be used as an unguided BDR mutation or as a
#' consensus-guided BDR mutation. When \code{consensus_score} is provided and
#' \code{consensus_lambda > 0}, the cross-island consensus information is used
#' to guide the locations involved in birth, death, and relocation moves.
#'
#' @details
#' For a chromosome
#' \deqn{C = (m, \boldsymbol{s}, \boldsymbol{\tau}, N+1)',}
#' where \eqn{m} is the number of changepoints,
#' \eqn{\boldsymbol{s}} contains optional model-order parameters, and
#' \eqn{\boldsymbol{\tau}} contains the ordered changepoint locations, the BDR
#' mutation applies exactly one feasible structural operation.
#'
#' A birth move adds one changepoint at a feasible location satisfying the
#' boundary and minimum-spacing constraints. A death move removes one existing
#' changepoint. A relocation move selects one existing changepoint and moves it
#' to another feasible location while keeping the number of changepoints fixed.
#'
#' When all three operations are feasible, each is selected with probability
#' one third. When only one or two operations are feasible, the probability is
#' redistributed equally among the available operations.
#'
#' If \code{consensus_score = NULL} or \code{consensus_lambda = 0}, the
#' operator performs unguided BDR mutation. Otherwise, consensus information
#' guides the mutation. Birth locations with stronger consensus support receive
#' larger sampling weights. During a death move, changepoints with weaker
#' consensus support receive larger deletion weights. During relocation,
#' weakly supported changepoints are more likely to be selected as the source,
#' while strongly supported feasible locations are more likely to be selected
#' as the destination.
#'
#' The model-order parameters in \eqn{\boldsymbol{s}}, when present, are not
#' modified by this mutation operator.
#'
#' @param child A vector or one-column matrix containing the chromosome to be
#' mutated.
#' @param prange The default value is \code{NULL} for changepoint detection only
#' task. If model order selection and changepoint detection are performed
#' simultaneously, \code{prange} should be a list specifying the allowable
#' ranges of the model-order parameters.
#' @param minDist The minimum length between two adjacent changepoints.
#' @param pchangepoint The probability that a changepoint has occurred. This
#' argument is retained for compatibility with the mutation-operator interface.
#' @param lmax The maximum possible length of the chromosome representation.
#' @param mmax The maximum number of changepoints allowed in the time series.
#' @param N The sample size of the time series.
#' @param consensus_score An optional numerical vector of length \code{N}
#' containing the consensus support score for each candidate location.
#' The default value is \code{NULL}, corresponding to unguided BDR mutation.
#' @param consensus_lambda A numerical value controlling the strength of
#' consensus guidance. A value of zero gives unguided BDR mutation. Larger
#' values give greater weight to the consensus score. The default is \code{0}.
#'
#' @return A one-column integer matrix containing the mutated chromosome.
#'
#' @references
#' Li, M. (2026). Structured and Consensus-Guided Island Model Genetic
#' Algorithm for Multiple Changepoint Detection. Manuscript submitted
#' for publication.
#'
#' @seealso
#' \code{\link{cptgaisl}},
#' \code{\link{cptgascisl}}
#'
#' @export
mutation_birth_death_relocate <- function(
    child,
    prange = NULL,
    minDist,
    pchangepoint,
    lmax,
    mmax,
    N,
    consensus_score = NULL,
    consensus_lambda = 0
) {
  
  child <- as.integer(child)
  
  plen <- length(prange)
  K <- child[1L]
  tau_start <- plen + 2L
  K_max <- min(mmax, lmax - plen - 2L)
  
  if (plen > 0L) {
    hyper <- child[2L:(plen + 1L)]
  } else {
    hyper <- integer(0)
  }
  
  if (K > 0L) {
    tau <- sort(child[tau_start:(tau_start + K - 1L)])
  } else {
    tau <- integer(0)
  }
  
  # birth operator
  birth_candidates <- integer(0)
  
  if (K < K_max) {
    if (K == 0L) {
      # current no cpt
      lower <- 1L + minDist
      upper <- N - minDist
      if (lower <= upper) {birth_candidates <- seq.int(lower, upper)}
    } else {
      lower <- c(1L + minDist, tau + minDist)
      upper <- c(tau - minDist, N - minDist)
      valid_intervals <- which(lower <= upper) # remove invalid intervals caused by two close cpts
      if (length(valid_intervals) > 0L) {
        birth_candidates <- unlist(lapply(valid_intervals, function(j) seq.int(lower[j], upper[j])),use.names = FALSE)
      }
    }
  }
  
  # globally relocate operator
  relocate_candidates <- vector("list", K)
  
  if (K > 0L) {
    for (j in seq_len(K)) {
      # build relocate neighborhood globally
      lower <- if (j == 1L) 1L + minDist else tau[j - 1L] + minDist
      upper <- if (j == K) N - minDist else tau[j + 1L] - minDist
      if (lower <= upper) {
        candidates_j <- seq.int(lower, upper)
        candidates_j <- candidates_j[candidates_j != tau[j]]
        if (length(candidates_j) > 0L) {
          relocate_candidates[[j]] <- candidates_j
        }
      }
    }
  }
  relocatable <- which(lengths(relocate_candidates) > 0L)
  
  # decide birth/death/relocate which one to operate equal prob
  
  possible_operations <- character(0)
  
  if (length(birth_candidates) > 0L) possible_operations <- c(possible_operations, "birth")
  if (K > 0L) possible_operations <- c(possible_operations, "death")
  if (length(relocatable) > 0L) possible_operations <- c(possible_operations, "relocate")
  
  if (length(possible_operations) == 0L) {
    return(matrix(child, nrow = lmax, ncol = 1L))
  }
  
  operation <- sample(possible_operations, 1L)
  
  if (operation == "birth") {
    if (is.null(consensus_score) || consensus_lambda <= 0) {
      new_tau <- birth_candidates[sample.int(length(birth_candidates), 1L)]
    } else {
      birth_weights <- (1 - consensus_lambda) + consensus_lambda * consensus_score[birth_candidates]
      new_tau <- birth_candidates[sample.int(length(birth_candidates), 1L, prob = birth_weights)]
    }
    tau <- sort(c(tau, new_tau))
  }
  
  if (operation == "death") {
    if (is.null(consensus_score) || consensus_lambda <= 0) {
      remove_j <- sample.int(K, 1L)
    } else {
      death_weights <- (1 - consensus_lambda) + consensus_lambda * (1 - consensus_score[tau])
      remove_j <- sample.int(K, 1L, prob = death_weights)
    }
    
    tau <- tau[-remove_j]
  }
  
  if (operation == "relocate") {
    if (is.null(consensus_score) || consensus_lambda <= 0) {
      relocate_j <- relocatable[sample.int(length(relocatable), 1L)]
    } else {
      source_weights <- (1 - consensus_lambda) + consensus_lambda * (1 - consensus_score[tau[relocatable]])
      relocate_j <- relocatable[sample.int(length(relocatable), 1L, prob = source_weights)]
    }
    candidates_j <- relocate_candidates[[relocate_j]]
    if (is.null(consensus_score) || consensus_lambda <= 0) {
      tau[relocate_j] <- candidates_j[sample.int(length(candidates_j), 1L)]
    } else {
      destination_weights <- (1 - consensus_lambda) + consensus_lambda * consensus_score[candidates_j]
      tau[relocate_j] <- candidates_j[sample.int(length(candidates_j), 1L, prob = destination_weights)]
    }
    tau <- sort(tau)
  }
  
  # Reconstruct chromosome
  # Birth:    K -> K+1
  # Death:    K -> K-1
  # Relocate: K -> K
  
  K <- length(tau)
  
  childMut <- matrix(0L, nrow = lmax, ncol = 1L) # return an matrix object to match
  childMut[1L, 1L] <- K
  
  if (plen > 0L) {
    childMut[2L:(plen + 1L), 1L] <- hyper # BDR NOT change hyper-parameters
  }
  
  if (K > 0L) {
    childMut[tau_start:(tau_start + K - 1L), 1L] <- tau
  }
  
  childMut[tau_start + K, 1L] <- N + 1L
  
  childMut
  
}