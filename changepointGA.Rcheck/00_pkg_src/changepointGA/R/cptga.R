#' Genetic algorithm
#'
#' Perform the genetic algorithm for changepoint detection.
#' This involves the minimization of an objective function using a genetic algorithm (GA).
#' The algorithm can be run sequentially or with explicit parallelization.
#'
#' For any pre-specified time series model with a specified set of changepoint locations,
#' model fit is evaluated using a fitness function \eqn{Q(\boldsymbol{\theta})},
#' where \eqn{\boldsymbol{\theta}=(\boldsymbol{s},\boldsymbol{\tau},\boldsymbol{\beta})'}
#' denotes the full parameter vector. Here, \eqn{\boldsymbol{s}} denotes the set
#' of model hyperparameters, potentially including the AR or MA orders, the degree
#' of ARIMA differencing, or the periodic AR order for PAR models. The vector
#' \eqn{\boldsymbol{\tau}=\{\tau_{1}, \ldots, \tau_{m}\}} specifies the
#' changepoint locations, with the number of changepoints \eqn{m} inferred as
#' part of the estimation. Each individual chromosome representation is specified as a vector,
#' \deqn{C = (m, \boldsymbol{s}, \boldsymbol{\tau}, N+1)',}
#' where \eqn{m} represents the number of changepoints and is also equivalent to
#' the length of vector \eqn{\boldsymbol{\tau}} containing all the changepoint
#' locations. \eqn{\boldsymbol{s}} contains the integer-valued parameter for time
#' series model structure specification, such as AR, MA, or PAR orders.
#' If no model order selection is desired, then \eqn{\boldsymbol{s}} is omitted
#' and the GA detects changepoints only. The changepoint locations in
#' \eqn{\boldsymbol{\tau}} are encoded as interger values between 2 and \eqn{N},
#' allowing the length of \eqn{\boldsymbol{\tau}} to vary dynamically with \eqn{m}.
#' The value \eqn{N+1} is appended at the end of \eqn{C} to serve as a delimiter
#' marking the boundary of the chromosome.
#'
#' @param ObjFunc The fitness function to be minimized. Users can specify any R or Rcpp
#' function as the fitness function, setting the input as the potential solution to
#' the optimization problem and returning a numerical value as the output/fitness.
#' Depending on the user-specified chromosome representation, the optimization task
#' can be changepoint detection only or changepoint detection plus model order selection,
#' which can be specified via the \code{option} parameter. When
#' \code{option="both"}, the list \code{prange} must be specified to give the range
#' of model orders.
#' @param N The sample size of the time series.
#' @param prange The default value is \code{NULL} for changepoint detection only task.
#' If both model order selection and changepoint detection are required, the \code{prange}
#' argument should be provided as a list. Each element in this list must specify the range
#' for a corresponding model order parameter, and the length of the list object of
#' \code{prange} must match the number of order parameters to be estimated.
#' @param popSize An integer represents the number of individuals/chromosomes in each population.
#' @param pcrossover The probability that the crossover operator will apply to
#' the two selected parents' chromosomes to produce the offspring. The typical
#' value is close to 1, with the default setting in this package being 0.95.
#' @param pmutation The probability that the mutation operator applies on one
#' individual chromosome. Similar to the natural mutation process, new genetic
#' information is introduced to the offspring chromosome with a relatively small
#' probability (close to 0), with a default value of 0.3.
#' @param pchangepoint The probability that a changepoint has occurred with a default value of 0.01. User
#' could change this probability based on domain knowledge and the time series
#' length. This probability is used during population initialization and in the
#' creation of new chromosomes by the mutation operator. By default, the mutation
#' operator function generates a new individual as the mutated offspring.
#' @param minDist The minimum length between two adjacent changepoints. Default value equals to one.
#' @param mmax The maximum number of changepoints allowed in the time series data
#' corresponds to the maximum possible length of \eqn{\boldsymbol{\tau}}.
#' For a time series of length 1000 and we only want to detect the changepoint
#' (\code{option="cp"}), the default value is 499. The suggested value should be
#' based on the length of the time series. For instance, if a time series has
#' length \eqn{N}, the recommended \code{mmax} should be \eqn{N/2-1}. It is suggested to
#' add the number of model hyperparameters if both changepoint detection and
#' model order selection tasks are of-interested simultaneously
#' (\code{option="both"}).
#' @param lmax The maximum possible length of the chromosome representation.
#' For a time series of length 1000 and we only want to detect the changepoint
#' (\code{option="cp"}), the default value is 501. The suggested value should be
#' based on the length of the time series. For instance, if a time series has
#' length \eqn{N}, the recommended \code{lmax} should be \eqn{2+N/2-1}. It is
#' suggested to add the number of model hyperparameters if both changepoint detection and
#' model order selection tasks are of-interested simultaneously (\code{option="both"}).
#' @param maxgen The maximum number of generations the GA can run before the
#' search is forcibly terminated.
#' @param maxconv If the overall best fitted value doesn't change after
#' \code{maxconv} consecutive migrations, the GA algorithm stops.
#' @param option A character string controls the optimization task. \code{"cp"}
#' indicates the task is changepoint detection only. \code{"both"} indicates the
#' task will include both changepoint detection and model order selection.
#' @param monitoring A logical value with either \code{TRUE} or \code{FALSE},
#' indicating whether to print out summarized results (current best fitness
#' function value and its corresponding $C$) for each generation from the GA.
#' @param parallel A logical value with either \code{TRUE} or \code{FALSE}, indicating
#' whether to use multiple cores for parallel computation of the fitness function
#' values for individuals after population initialization.
#' @param nCore An integer with the default value of \code{NULL}. It
#' represents the number of cores used in parallel computing. The value of
#' \code{nCore} must be less than the number of physical cores available.
#' @param tol An numerical value with the default value of \code{1e-05}. The
#' tolerance level that helps evaluate whether the two iterations have the same
#' fitness value, which aids in determining GA termination.
#' @param seed An integer with the default value of \code{NULL}. An single
#' integer allows function produce reproducible results.
#' @param popInitialize A function or sourced function name character string.
#' It should be designed for initializing a population. The default population
#' initialization is random initialization with some imposed constraints. See
#' \code{\link{random_population}} for example. The function returned object is
#' a matrix, \code{pop}. The users can specified their own \code{population}
#' function. It could also be a matrix object, which contain the user specified
#' chromosome. By default, each column represents one individual chromosome.
#' See \code{\link{random_population}} for details.
#' @param suggestions A list object. Default value is \code{NULL}. Each element
#' includes better or more reasonable guessed changepoints locations, which will
#' resulting in one chromosome. If the number of provided suggestions equals to
#' \code{popSize}, all the suggestions will be used as population. If the number
#' of provided suggestions is less than \code{popSize}, the function from
#' \code{popInitialize} will generate the remaining individual chromosomes.
#' The number of provided suggestions cannot be greater than \code{popSize}.
#' Having better \code{suggestions} can help GA converges faster.
#' @param selection A function or sourced function name character string. This
#' GA operator can help select \code{mom} and \code{dad} from current generation
#' population, where \code{dad} is set to have better fit (smaller fitness
#' function values). The default for selection uses the linear rank selection
#' method. See \code{\link{selection_linearrank}} for example. The function
#' returned object is a list contain the chromosomes for \code{mom} and
#' \code{dad}.
#' @param crossover A function or sourced function name character string. This
#' GA operator can apply crossover to the chosen parents to produce child for
#' next generation with specified probability. The default for crossover uses
#' the uniform crossover method. See \code{\link{uniformcrossover}} for details
#' in the default crossover operator. The function returned object is a vector
#' contain the chromosomes for \code{child}.
#' @param mutation A function or sourced function name character string. This
#' GA operator can apply mutation to the produced \code{child} with the
#' specified probability \code{pmutation}. See \code{\link{mutation}} for
#' details in the default mutation operator. The function returned object
#' is a vector contain \code{child} chromosome representation.
#' @param ... additional arguments that will be passed to the fitness function.
#' @return Return an object class \code{cptga-class}. See
#' \code{\link{cptga-class}} for a more detailed description.
#' @import stats
#' @import Rcpp
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom methods new
#' @importFrom utils str
#' @useDynLib changepointGA
#' @export
#' @examples
#' \donttest{
#'
#' N <- 1000
#' XMatT <- matrix(1, nrow = N, ncol = 1)
#' Xt <- ts.sim(
#'   beta = 0.5, XMat = XMatT, sigma = 1, phi = 0.5, theta = NULL,
#'   Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
#' )
#'
#' ## Multiple changepoint detection without model order selection
#'
#' # without suggestions
#' GA.res <- cptga(ObjFunc = ARIMA.BIC, N = N, XMat = XMatT, Xt = Xt)
#' summary(GA.res)
#' plot(GA.res, data = Xt)
#'
#' # with suggestions
#' suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
#' GA.res <- cptga(ObjFunc = ARIMA.BIC, N = N, suggestions = suggestions, XMat = XMatT, Xt = Xt)
#' summary(GA.res)
#' plot(GA.res, data = Xt)
#'
#'
#' ## Multiple changepoint detection with model order selection
#'
#' prange <- list(ar = c(0, 3), ma = c(0, 3))
#'
#' # without suggestions
#' GA.res <- cptga(
#'   ObjFunc = ARIMA.BIC.Order, N = N, prange = prange,
#'   option = "both", XMat = XMatT, Xt = Xt
#' )
#' summary(GA.res)
#' plot(GA.res, data = Xt)
#'
#' # with suggestions
#' suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
#' GA.res <- cptga(
#'   ObjFunc = ARIMA.BIC.Order, N = N, prange = prange,
#'   suggestions = suggestions, option = "both", XMat = XMatT, Xt = Xt
#' )
#' summary(GA.res)
#' plot(GA.res, data = Xt)
#' }
cptga <- function(ObjFunc,
                  N,
                  prange = NULL,
                  popSize = 200,
                  pcrossover = 0.95,
                  pmutation = 0.3,
                  pchangepoint = 0.01,
                  minDist = 1,
                  mmax = NULL,
                  lmax = NULL,
                  maxgen = 50000,
                  maxconv = 5000,
                  option = "cp",
                  monitoring = FALSE,
                  parallel = FALSE,
                  nCore = NULL,
                  tol = 1e-05,
                  seed = NULL,
                  popInitialize = "random_population",
                  suggestions = NULL,
                  selection = "selection_linearrank",
                  crossover = "uniformcrossover",
                  mutation = "mutation",
                  ...) {
  call <- match.call()
  dots <- list(...)

  if (missing(N)) {
    stop("The sample size must be provided")
  }

  if (is.null(mmax)) {
    mmax <- floor(N / 2 - 1)
  }
  if (is.null(lmax)) {
    lmax <- floor(2 + N / 2 - 1)
  }

  if (missing(ObjFunc)) {
    stop("A fitness function must be provided")
  }
  if (!is.function(ObjFunc)) {
    stop("A fitness function must be provided")
  }

  i <- NULL # add for global variables declare
  plen <- length(prange)

  if (pcrossover < 0 | pcrossover > 1) {
    stop("Probability of crossover must be between 0 and 1.")
  }
  if (pmutation < 0 | pmutation > 1) {
    stop("Probability of mutation must be between 0 and 1.")
  }
  if (pchangepoint < 0 | pchangepoint > 1) {
    stop("Probability of changepoint must be between 0 and 1.")
  }
  if (minDist >= N | minDist < 1) {
    stop("Minimum number of locations between two changepoints invalid.")
  }
  if (lmax < mmax + 2) {
    stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).")
  }
  if (option == "both" & plen == 0) {
    stop("Opt for changepoint and order search, prange must be provided.")
  }
  if (option == "cp" & plen != 0) {
    stop("Opt for changepoint search, prange must be NULL.")
  }

  # set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.function(popInitialize)) {
    popInitialize <- get(popInitialize)
  }
  if (!is.function(selection)) {
    selection <- get(selection)
  }
  if (!is.function(crossover)) {
    crossover <- get(crossover)
  }
  if (!is.function(mutation)) {
    mutation <- get(mutation)
  }

  object <- new("cptga",
    call = call,
    N = N,
    prange = prange,
    popSize = popSize,
    pcrossover = pcrossover,
    pmutation = pmutation,
    pchangepoint = pchangepoint,
    minDist = minDist,
    mmax = mmax,
    lmax = lmax,
    maxgen = maxgen,
    maxconv = maxconv,
    option = option,
    monitoring = monitoring,
    parallel = parallel,
    nCore = nCore,
    tol = tol,
    seed = seed,
    suggestions = suggestions,
    population = matrix(),
    fitness = vector(),
    overbestchrom = vector(),
    overbestfit = double(),
    bestfit = vector(),
    count = integer(),
    convg = integer()
  )

  ###### step 1: Initialize population
  if (!is.null(suggestions)) {
    if (!is.list(suggestions)) {
      stop("Suggested changepoints must be a List.")
    }
    # 1. with suggestions
    if (any(sapply(suggestions, function(x) any(is.na(x))))) {
      stop("NA provided in suggestions list")
    }
    n.suggs <- length(suggestions)
    suggestions.mat <- matrix(0L, nrow = lmax, ncol = n.suggs)
    for (i in seq_along(suggestions)) {
      idx <- suggestions[[i]]
      if (any(idx <= 1 | idx > N)) {
        stop(paste0("\n No. ", i, "th suggestion is invalid."))
      }
      if (object@option == "cp") {
        if (is.null(idx)) {
          suggestions.mat[1:2, i] <- c(0, N + 1)
        } else {
          idx <- idx[idx > 1 & idx <= N] # filter valid indices
          suggestions.mat[1:(length(idx) + 2), i] <- c(length(idx), idx, N + 1)
        }
      } else if (object@option == "both") {
        hyper.param <- sapply(prange, function(x) if (length(x) > 0) sample(min(x):max(x), 1) else NA)
        if (is.null(idx)) {
          suggestions.mat[1:(2 + plen), i] <- c(0, hyper.param, N + 1)
        } else {
          idx <- idx[idx > 1 & idx <= N] # filter valid indices
          suggestions.mat[1:(length(idx) + plen + 2), i] <- c(length(idx), hyper.param, idx, N + 1)
        }
      }
    }
    RemainpopSize <- popSize - n.suggs
    if (RemainpopSize > 0) {
      # 1.1 partial from suggestions
      init_args <- c(
        dots,
        list(
          popSize = RemainpopSize, prange = prange, N = N,
          minDist = minDist, pchangepoint = pchangepoint,
          mmax = mmax, lmax = lmax
        )
      )
      init_args <- .filter_args(init_args, popInitialize)
      pop <- do.call(popInitialize, init_args)
      pop <- cbind(pop, suggestions.mat)
    } else if (RemainpopSize == 0) {
      # 1.2 completely from suggestions
      pop <- suggestions.mat
    } else {
      # 1.3 some errors
      stop("Number of suggested chromosome must be less than specified population size.")
    }
  } else {
    # 2. completely generated by popInitialize
    init_args <- c(
      dots,
      list(
        popSize = popSize, prange = prange, N = N,
        minDist = minDist, pchangepoint = pchangepoint,
        mmax = mmax, lmax = lmax
      )
    )
    init_args <- .filter_args(init_args, popInitialize)
    pop <- do.call(popInitialize, init_args)
  }

  obj_shared <- .filter_args(dots, ObjFunc)
  obj_formals <- setdiff(names(formals(ObjFunc)), "...")
  first_param <- obj_formals[1]
  has_plen <- "plen" %in% obj_formals

  ## evaluate the fitness (Parallel or NOT)
  if (parallel) {
    nAvaCore <- detectCores()
    if (is.null(nCore)) {
      stop(paste0("Missing number of computing cores (", nAvaCore, " cores available)."))
    }
    registerDoParallel(cores = nCore)
    popFit <- foreach::foreach(j = 1:popSize, .combine = "c") %dopar% {
      chrom <- pop[1:(pop[1, j] + plen + 2), j]
      call_list <- c(
        setNames(list(chrom), first_param),
        if (has_plen) list(plen = plen) else list(),
        obj_shared
      )
      do.call(ObjFunc, call_list)
    }
  } else {
    popFit <- rep(NA, popSize)
    for (j in 1:popSize) {
      chrom <- pop[1:(pop[1, j] + plen + 2), j]
      call_list <- c(
        setNames(list(chrom), first_param),
        if (has_plen) list(plen = plen) else list(),
        obj_shared
      )
      popFit[j] <- do.call(ObjFunc, call_list)
    }
  }

  object@population <- pop
  object@fitness <- popFit

  count <- 0
  bestfit <- rep(NA, maxgen)
  bestchrom <- matrix(0, nrow = lmax, ncol = maxgen)
  repeat{
    # indicator for c("crossover", "mutation")
    #     flag[1]=1 indicating no cross-over
    #     flag[2]=1 indicating no mutation
    flag <- rep(0, 2)

    ##### step 2: parents selection
    sel_args <- .filter_args(dots, selection)
    parents <- do.call(selection, c(list(pop, popFit), sel_args))
    dad <- parents$dad
    mom <- parents$mom

    ##### step 3: crossover
    a1 <- runif(1)
    if (a1 <= pcrossover) {
      cross_args <- .filter_args(dots, crossover)
      child <- do.call(crossover, c(
        list(
          mom = mom, dad = dad,
          prange = prange, minDist = minDist,
          lmax = lmax, N = N
        ),
        cross_args
      ))
    } else {
      child <- dad
      flag[1] <- 1
    }

    ##### step 4: mutation
    a2 <- runif(1)
    if (a2 <= pmutation) {
      mut_args <- .filter_args(dots, mutation)
      child <- do.call(mutation, c(
        list(
          child, prange, minDist,
          pchangepoint, lmax, mmax, N
        ),
        mut_args
      ))
    } else {
      flag[2] <- 1
    }

    ##### step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum <- flag[1] + flag[2]

    if (flagsum < 2) {
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      child_chrom <- child[1:(child[1] + plen + 2)]
      call_list <- c(
        setNames(list(child_chrom), first_param),
        if (has_plen) list(plen = plen) else list(),
        obj_shared
      )
      fitChild <- do.call(ObjFunc, call_list)
      leastfit <- max(popFit) # with largest fitness value

      if (fitChild < leastfit) {
        # indicating child is better than the worst one and replace
        pp <- which.max(popFit)
        pop[, pp] <- child
        popFit[pp] <- fitChild
      }
    }

    object@population <- pop
    object@fitness <- popFit

    # best popFit with minimized fitness value
    count <- count + 1
    genbest <- which.min(popFit)
    bestfit[count] <- popFit[genbest]
    bestchrom[, count] <- pop[, genbest]

    object@bestfit <- bestfit

    # check convergence once count >= maxconv
    # if after maxconv consecutive generations, the overall best does not change, then stop
    if (count >= maxconv) {
      tmpbestfit <- bestfit[(count - maxconv + 1):count]
      decision <- checkConv(tmpbestfit, maxconv, tol)
      if (monitoring) {
        cat("\n My decision:", decision)
      }
      if (decision == 1) {
        overbestfit <- bestfit[count]
        overbestchrom <- bestchrom[, count]
        overbestchrom <- overbestchrom[1:(overbestchrom[1] + plen + 2)]
        if (monitoring) {
          cat("\n==== No.", count, "Generation ====")
          cat("\n overall bestfit =", overbestfit)
          cat("\n overall bestchrom =", overbestchrom, "\n")
        }
        pp <- which(is.na(bestfit))[1] - 1
        bestfit <- bestfit[1:pp]
        bestchrom <- bestchrom[, 1:pp]

        object@count <- count
        object@convg <- 0
        object@bestfit <- bestfit

        break
      }
    }

    overbestfit <- bestfit[count]
    overbestchrom <- bestchrom[, count]
    object@overbestfit <- bestfit[count]
    object@overbestchrom <- overbestchrom[1:(overbestchrom[1] + plen + 2)]
    if (monitoring) {
      cat("\n==== No.", count, "Generation ====")
      cat("\n overall bestfit =", overbestfit)
      cat("\n overall bestchrom =", overbestchrom, "\n")
    }

    if (count >= maxgen) {
      object@count <- count
      object@convg <- 1
      break
    }
  }

  return(object)
}
