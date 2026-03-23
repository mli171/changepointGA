#' Island model based genetic algorithm
#'
#' Perform the modified island-based genetic algorithm (IslandGA) for multiple changepoint detection.
#' Minimization of an objective function using genetic algorithm (GA).
#' The algorithm can be run sequentially or in explicit parallelisation.
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
#' @param popSize An integer representing the total number of individuals in each generation, which
#' equal to the number of islands multiplied by the size of each island (i.e.,
#' \code{popSize = numIslands × Islandsize}).
#' @param numIslands An integer representing the number of islands (sub-populations).
#' @param pcrossover The probability that the crossover operator will apply to
#' the two selected parents' chromosomes to produce the offspring. The typical
#' value is close to 1, with the default setting in this package being 0.95.
#' @param pmutation The probability that the mutation operator applies on one
#' individual chromosome. Similar to the natural mutation process, new genetic
#' information is introduced to the offspring chromosome with a relatively small
#' probability (close to 0), with a default value of 0.3.
#' @param pchangepoint The probability that a changepoint has occurred. User
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
#' length \eqn{N}, the recommended \code{lmax} should be \eqn{2+N/2-1}. It is suggested to
#' add the number of model hyperparameters if both changepoint detection and
#' model order selection tasks are of-interested simultaneously (\code{option="both"}).
#' @param maxMig An integer indicates the maximum number of migrations. After
#' conducting \code{maxMig} migrations, the island model GA stops.
#' @param maxgen An integer indicates the maximum number of generations that each
#' island (subpopulation) undergoes before migration. It also determines the
#' requency of migration. The migration will occur after \code{maxgen} generations
#' for each island (sub-population).
#' @param maxconv An integer value is also used for algorithm termination. If
#' the overall best-fitted value remains unchanged after \code{maxconv}
#' consecutive migrations, the island model GA will terminate.
#' @param option A character string controls the optimization task. \code{"cp"}
#' indicates the task is changepoint detection only. \code{"both"} indicates the
#' task will include both changepoint detection and model order selection.
#' @param monitoring A logical value with either \code{TRUE} or \code{FALSE},
#' indicating whether to print out summarized results (current best fitness
#' function value and its corresponding $C$) for each generation from the GA.
#' @param parallel A logical value with \code{TRUE} or \code{FALSE}, indicating
#' whether to use multiple cores for parallel computation. Different from the
#' basic GA, in this approach, each island evolves independently. New population
#' generation—including parent selection, crossover, and mutation—is processed
#' independently and in parallel for each island, further saving computation
#' time. Once the new population generation is complete, the migration operator
#' is performed sequentially.
#' @param nCore An integer indicates the number of cores utilized in parallel
#' computing. For the island model GA, the number of cores used in parallel is
#' set to the \code{numIslands} for convenience.
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
#' method. See \code{\link{selection_linear_rank}} for example. The function
#' returned object is a list contain the chromosomes for \code{mom} and
#' \code{dad}.
#' @param crossover A function or sourced function name character string. This
#' GA operator can apply crossover to the chosen parents to produce child for
#' next generation with specified probability. The default for crossover uses
#' the uniform crossover method. See \code{\link{uniform_crossover}} for details
#' in the default crossover operator. The function returned object is a vector
#' contain the chromosomes for \code{child}.
#' @param mutation A function or sourced function name character string. This
#' GA operator can apply mutation to the produced \code{child} with the
#' specified probability \code{pmutation}. See \code{\link{mutation}} for
#' details in the default mutation operator. The function returned object
#' is a vector contain \code{child} chromosome representation.
#' @param ... additional arguments that will be passed to the fitness function.
#' @return Return an object class \code{cptgaisl-class}. See
#' \code{\link{cptgaisl-class}} for a more detailed description.
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
#' Xt <- ts_sim(
#'   beta = 0.5, XMat = XMatT, sigma = 1, phi = 0.5, theta = NULL,
#'   Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
#' )
#'
#' ## Multiple changepoint detection without model order selection
#'
#' # without suggestions
#' GAISLres <- cptgaisl(ObjFunc = arima_bic, N = N, XMat = XMatT, Xt = Xt)
#' summary(GAISLres)
#' plot(GAISLres, data = Xt)
#'
#' # with suggestions
#' suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
#' GAISLres <- cptgaisl(ObjFunc = arima_bic, N = N, suggestions = suggestions, XMat = XMatT, Xt = Xt)
#' summary(GAISLres)
#' plot(GAISLres, data = Xt)
#'
#'
#' ## Multiple changepoint detection with model order selection
#'
#' prange <- list(ar = c(0, 3), ma = c(0, 3))
#'
#' # without suggestions
#' GAISLres <- cptgaisl(
#'   ObjFunc = arima_bic_order, N = N, prange = prange,
#'   option = "both", XMat = XMatT, Xt = Xt
#' )
#' summary(GAISLres)
#' plot(GAISLres, data = Xt)
#'
#' # with suggestions
#' suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
#' GAISLres <- cptgaisl(
#'   ObjFunc = arima_bic_order, N = N, prange = prange,
#'   suggestions = suggestions, option = "both", XMat = XMatT, Xt = Xt
#' )
#' summary(GAISLres)
#' plot(GAISLres, data = Xt)
#' }
cptgaisl <- function(ObjFunc,
                     N,
                     prange = NULL,
                     popSize = 200,
                     numIslands = 5,
                     pcrossover = 0.95,
                     pmutation = 0.3,
                     pchangepoint = 0.01,
                     minDist = 1,
                     mmax = NULL,
                     lmax = NULL,
                     maxMig = 1000,
                     maxgen = 50,
                     maxconv = 100,
                     option = "cp",
                     monitoring = FALSE,
                     parallel = FALSE,
                     nCore = NULL,
                     tol = 1e-05,
                     seed = NULL,
                     popInitialize = "random_population",
                     suggestions = NULL,
                     selection = "selection_linear_rank",
                     crossover = "uniform_crossover",
                     mutation = "mutation",
                     ...) {
  call <- match.call()
  dots <- list(...)

  if (missing(N)) {
    stop("The sample size must be provided")
  }
  if (missing(ObjFunc) || !is.function(ObjFunc)) {
    stop("`ObjFunc` must be a function.")
  }
  
  if (is.null(mmax)) {
    mmax <- floor(N / 2 - 1)
  }
  if (is.null(lmax)) {
    lmax <- floor(2 + N / 2 - 1)
  }
  
  plen <- length(prange)

  if (pcrossover < 0 || pcrossover > 1) {
    stop("Probability of crossover must be between 0 and 1.")
  }
  if (pmutation < 0 || pmutation > 1) {
    stop("Probability of mutation must be between 0 and 1.")
  }
  if (pchangepoint < 0 || pchangepoint > 1) {
    stop("Probability of changepoint must be between 0 and 1.")
  }
  if (minDist >= N || minDist < 1) {
    stop("Minimum number of locations between two changepoints invalid.")
  }
  if (lmax < mmax + 2) {
    stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).")
  }
  if (option == "both" && plen == 0) {
    stop("Opt for changepoint and order search, prange must be provided.")
  }
  if (option == "cp" && plen != 0) {
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

  object <- new("cptgaisl",
    call = call,
    N = N,
    prange = prange,
    popSize = popSize,
    numIslands = numIslands,
    Islandsize = double(),
    pcrossover = pcrossover,
    pmutation = pmutation,
    pchangepoint = pchangepoint,
    minDist = minDist,
    mmax = mmax,
    lmax = lmax,
    maxMig = maxMig,
    maxgen = maxgen,
    maxconv = maxconv,
    option = option,
    monitoring = monitoring,
    parallel = parallel,
    nCore = nCore,
    tol = tol,
    seed = seed,
    suggestions = suggestions,
    Island = array(),
    IslandFit = matrix(),
    overbestchrom = vector(),
    overbestfit = double(),
    bestfit = vector(),
    countMig = integer(),
    count = integer(),
    convg = integer()
  )

  ###### step 1: generate initial population for each island
  if (!is.null(suggestions)) {
    if (!is.list(suggestions)) {
      stop("`suggestions` must be a list.")
    }
    if (any(vapply(suggestions, function(x) any(is.na(x)), logical(1)))) {
      stop("NA values are not allowed in `suggestions`.")
    }
    
    n_suggs <- length(suggestions)
    suggestions_mat <- matrix(0L, nrow = lmax, ncol = n_suggs)
    
    for (i in seq_along(suggestions)) {
      idx <- suggestions[[i]]
      
      if (!is.null(idx) && any(idx <= 1 | idx > N)) {
        stop("Suggestion ", i, " is invalid.")
      }
      
      if (option == "cp") {
        if (is.null(idx)) {
          suggestions_mat[1:2, i] <- c(0L, N + 1L)
        } else {
          idx <- idx[idx > 1 & idx <= N]
          suggestions_mat[1:(length(idx) + 2L), i] <- c(length(idx), idx, N + 1L)
        }
      } else {
        hyper_param <- vapply(
          prange,
          function(x) sample(min(x):max(x), 1),
          integer(1)
        )
        
        if (is.null(idx)) {
          suggestions_mat[1:(plen + 2L), i] <- c(0L, hyper_param, N + 1L)
        } else {
          idx <- idx[idx > 1 & idx <= N]
          suggestions_mat[1:(length(idx) + plen + 2L), i] <- c(length(idx), hyper_param, idx, N + 1L)
        }
      }
    }
    
    remain_pop_size <- popSize - n_suggs
    
    if (remain_pop_size > 0) {
      init_args <- c(
        dots,
        list(
          popSize = remain_pop_size,
          prange = prange,
          N = N,
          minDist = minDist,
          pchangepoint = pchangepoint,
          mmax = mmax,
          lmax = lmax
        )
      )
      init_args <- .filter_args(init_args, popInitialize)
      pop <- do.call(popInitialize, init_args)
      pop <- cbind(pop, suggestions_mat)
    } else if (remain_pop_size == 0) {
      pop <- suggestions_mat
    } else {
      stop("The number of suggested chromosomes must not exceed `popSize`.")
    }
  } else {
    init_args <- c(
      dots,
      list(
        popSize = popSize,
        prange = prange,
        N = N,
        minDist = minDist,
        pchangepoint = pchangepoint,
        mmax = mmax,
        lmax = lmax
      )
    )
    init_args <- .filter_args(init_args, popInitialize)
    pop <- do.call(popInitialize, init_args)
  }

  Islandsize <- floor(popSize / numIslands)
  remainder <- popSize %% numIslands
  # Shuffle individual indices
  shuffled_idx <- sample.int(popSize)
  
  Island <- array(0L, dim = c(lmax, Islandsize, numIslands))
  for (k in seq_len(numIslands)) {
    idx_start <- (k - 1L) * Islandsize + 1L
    idx_end <- k * Islandsize
    idx <- shuffled_idx[idx_start:idx_end]
    Island[, , k] <- pop[, idx, drop = FALSE]
  }
  
  if (remainder > 0L) {
    warning(remainder, " individuals were unassigned. Consider adjusting `numIslands`.")
    popSize <- numIslands * Islandsize
    object@popSize <- popSize
  }

  object@numIslands <- numIslands
  object@Islandsize <- Islandsize
  
  obj_shared <- .filter_args(dots, ObjFunc)
  obj_formals <- setdiff(names(formals(ObjFunc)), "...")
  first_param <- obj_formals[1]
  has_plen <- "plen" %in% obj_formals

  eval_fitness <- function(chromosome) {
    chrom <- trim_chromosome(chromosome, plen)
    call_list <- c(
      setNames(list(chrom), first_param),
      if (has_plen) list(plen = plen) else list(),
      obj_shared
    )
    out <- do.call(ObjFunc, call_list)
    if (!is.finite(out)) Inf else out
  }
  
  get_island_best <- function(pop_mat, fit_vec) {
    best_idx <- which.min(fit_vec)
    list(
      fit = fit_vec[best_idx],
      chrom = pop_mat[, best_idx]
    )
  }
  
  ## evaluate the fitness (Parallel or NOT)
  if (parallel) {
    n_ava_core <- parallel::detectCores()
    if (is.null(nCore)) {
      nCore <- numIslands
    }
    if (nCore < 1 || nCore > n_ava_core) {
      stop("`nCore` must be between 1 and ", n_ava_core, ".")
    }
    
    doParallel::registerDoParallel(cores = nCore)
    on.exit(foreach::registerDoSEQ(), add = TRUE)
    
    IslandFit <- foreach::foreach(k = seq_len(numIslands), .combine = "cbind") %dopar% {
      vapply(seq_len(Islandsize), function(j) eval_fitness(Island[, j, k]), numeric(1))
    }
  } else {
    IslandFit <- matrix(0, nrow = Islandsize, ncol = numIslands)
    for (k in seq_len(numIslands)) {
      IslandFit[, k] <- vapply(seq_len(Islandsize), function(j) eval_fitness(Island[, j, k]), numeric(1))
    }
  }

  object@Island <- Island
  object@IslandFit <- IslandFit

  countMig <- 0L
  countMig <- 0L
  Bfit <- rep(Inf, numIslands)
  Bchrom <- matrix(0L, nrow = lmax, ncol = numIslands)
  bestfit <- rep(NA_real_, maxMig)
  bestchrom <- matrix(0L, nrow = lmax, ncol = maxMig)

  while (countMig < maxMig) {
    
    # step 2,3,4,5: select parents, crossover, mutation, new pop
    if (parallel) {
      resNewpop <- foreach::foreach(k = seq_len(numIslands)) %dopar% {
        do.call(new_population_Island, c(
          dots,
          list(
            ObjFunc = ObjFunc,
            prange = prange,
            selection = selection,
            crossover = crossover,
            mutation = mutation,
            pop = Island[, , k],
            fit = IslandFit[, k],
            minDist = minDist,
            lmax = lmax,
            mmax = mmax,
            pcrossover = pcrossover,
            pmutation = pmutation,
            pchangepoint = pchangepoint,
            maxgen = maxgen,
            N = N
          )
        ))
      }
      
      for (k in seq_len(numIslands)) {
        tmpfit <- resNewpop[[k]][1, ]
        tmppop <- resNewpop[[k]][-1, , drop = FALSE]
        tmpfit[!is.finite(tmpfit)] <- Inf
        
        Island[, , k] <- tmppop
        IslandFit[, k] <- tmpfit
        
        island_best <- get_island_best(tmppop, tmpfit)
        Bfit[k] <- island_best$fit
        Bchrom[, k] <- island_best$chrom
      }
    } else {
      for (k in seq_len(numIslands)) {
        resNewpop <- do.call(new_population_Island, c(
          dots,
          list(
            ObjFunc = ObjFunc,
            prange = prange,
            selection = selection,
            crossover = crossover,
            mutation = mutation,
            pop = Island[, , k],
            fit = IslandFit[, k],
            minDist = minDist,
            lmax = lmax,
            mmax = mmax,
            pcrossover = pcrossover,
            pmutation = pmutation,
            pchangepoint = pchangepoint,
            maxgen = maxgen,
            N = N
          )
        ))
        tmpfit <- resNewpop[1, ]
        tmppop <- resNewpop[-1, , drop = FALSE]
        tmpfit[!is.finite(tmpfit)] <- Inf
        
        Island[, , k] <- tmppop
        IslandFit[, k] <- tmpfit
        
        island_best <- get_island_best(tmppop, tmpfit)
        Bfit[k] <- island_best$fit
        Bchrom[, k] <- island_best$chrom
      }
    }
    
    # step 6: migration
    for (k in seq_len(numIslands)) {
      worst_idx <- which.max(IslandFit[, k])
      donor_pool <- setdiff(seq_len(numIslands), k)
      if (length(donor_pool) == 0L) {next}
      donor_island <- sample(donor_pool, 1)
      Island[, worst_idx, k] <- Bchrom[, donor_island]
      IslandFit[worst_idx, k] <- Bfit[donor_island]
    }
    ## recompute island bests after migration
    for (k in seq_len(numIslands)) {
      island_best <- get_island_best(Island[, , k], IslandFit[, k])
      Bfit[k] <- island_best$fit
      Bchrom[, k] <- island_best$chrom
    }
    
    # update best
    countMig <- countMig + 1L
    genbest <- which.min(Bfit)
    bestfit[countMig] <- Bfit[genbest]
    bestchrom[, countMig] <- Bchrom[, genbest]
    
    object@countMig <- countMig
    object@count <- countMig * maxgen
    object@Island <- Island
    object@IslandFit <- IslandFit
    object@bestfit <- bestfit[seq_len(countMig)]
    
    if (countMig >= maxconv) {
      tmpbest <- bestfit[(countMig - maxconv + 1L):countMig]
      decision <- check_conv(tmpbest, maxconv, tol)
      
      if (monitoring) {
        message("Decision: ", decision)
      }
      
      if (decision == 1L) {
        object@convg <- 0L
        object@overbestfit <- bestfit[countMig]
        object@overbestchrom <- trim_chromosome(bestchrom[, countMig], plen)
        return(object)
      }
    }
    
    if (monitoring) {
      current_bestfit <- bestfit[countMig]
      current_bestchrom <- trim_chromosome(bestchrom[, countMig], plen)
      
      message(
        "==== Migration ", countMig, " ====\n",
        "countMig = ", countMig, "\n",
        "overbestfit = ", current_bestfit, "\n",
        "overbestchrom = ", paste(current_bestchrom, collapse = " ")
      )
    }
  }
  
  object@convg <- 1L
  object@overbestfit <- bestfit[countMig]
  object@overbestchrom <- trim_chromosome(bestchrom[, countMig], plen)
  
  return(object)
  
}
