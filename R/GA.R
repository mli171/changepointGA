#' Genetic algorithm
#'
#' Perform the modified genetic algorithm for changepoint detection.
#' This involves the minimization of an objective function using a genetic algorithm (GA).
#' The algorithm can be run sequentially or with explicit parallelization.
#'
#' @param ObjFunc The fitness function to be minimized. Users can specify any R or Rcpp
#' function as the fitness function, setting the input as the potential solution to
#' the optimization problem and returning a numerical value as the output/fitness.
#' Depending on the user-specified chromosome representation, the optimization task
#' can be changepoint detection only or changepoint detection plus model order selection,
#' which can be specified via the \code{option} parameter in \code{\link{GA_param}}. When
#' \code{option="both"}, the list \code{p.range} must be specified to give the range
#' of model orders.
#' @param N The sample size of the time series.
#' @param GA_param A list contains the hyper-parameters for genetic algorithm.
#' The default values for these hyper-parameters are included in \code{.default.GA_param}.
#' See \code{\link{GA_param}} for the details.
#' @param GA_operators A list includes the functions for population initialization,
#' new individual selection, and genetic operator of crossover and mutation.
#' See \code{\link{operators}} for the details and default functions.
#' @param p.range Default is \code{NULL} for only changepoint detection. If
#' \code{p.range} is specified as a list object, which contains the range of
#' each model order parameters for order selection (integers). The number of
#' order parameters must be equal to the length of \code{p.range}.
#' @param ... additional arguments that will be passed to the fitness function.
#' @return Returns a list that has components:
#' \item{overbestfit}{The obtained minimum value of the objective function after
#' the final iteration.}
#' \item{overbestchrom}{The locations of the detected changepoints associating
#' with the \code{overbestfit} the after the final iteration.}
#' \item{bestfit}{The minimized fitness function values at each iteration.}
#' \item{bestchrom}{The detected changepoints at each iteration.}
#' \item{count}{The number of iterations undertaken by the genetic algorithm.}
#' \item{convg}{An integer code.
#'  \itemize{
#'    \item{0} indicates the algorithm successful completion.
#'    \item{1} indicates the the total number of generations exceeds the
#'             prespecified \code{maxgen} limit.
#'  }
#' }
#'
#' @import stats
#' @import Rcpp
#' @import foreach
#' @import doParallel
#' @import RcppArmadillo
#' @import parallel
#' @useDynLib changepointGA
#' @export
#' @examples
#' \donttest{
#' N = 1000
#' XMatT = matrix(1, nrow=N, ncol=1)
#' Xt = ts.sim(beta=0.5, XMat=XMatT, sigma=1, phi=0.5, theta=NULL,
#'             Delta=c(2, -2), CpLoc=c(250, 750), seed=1234)
#' TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=c(250, 750))
#'
#' GA.res = GA(ObjFunc=ARIMA.BIC, N=N, XMat=XMatT, Xt=Xt)
#' GA.res$overbestfit
#' GA.res$overbestchrom
#' }
GA = function(ObjFunc, N, GA_param=.default.GA_param, GA_operators=.default.operators, p.range=NULL, ...){

  call = match.call()

  i = NULL # add for global variables declare
  plen = length(p.range)

  popsize      = GA_param$popsize
  Pcrossover   = GA_param$Pcrossover
  Pmutation    = GA_param$Pmutation
  Pchangepoint = GA_param$Pchangepoint
  minDist      = GA_param$minDist
  mmax         = GA_param$mmax
  lmax         = GA_param$lmax
  maxgen       = GA_param$maxgen
  maxconv      = GA_param$maxconv
  option       = GA_param$option
  monitoring   = GA_param$monitoring
  parallel     = GA_param$parallel
  nCore        = GA_param$nCore
  tol          = GA_param$tol
  seed         = GA_param$seed

  if(missing(ObjFunc))
  { stop("A fitness function must be provided") }
  if(!is.function(ObjFunc))
  { stop("A fitness function must be provided") }
  if(Pcrossover < 0   | Pcrossover > 1)
  { stop("Probability of crossover must be between 0 and 1.") }
  if(Pmutation < 0    | Pmutation > 1)
  { stop("Probability of mutation must be between 0 and 1.") }
  if(Pchangepoint < 0 | Pchangepoint > 1)
  { stop("Probability of changepoint must be between 0 and 1.") }
  if(minDist >= N | minDist < 1)
  { stop("Minimum number of locations between two changepoints invalid.") }
  if(lmax < mmax + 2 )
  { stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).") }
  if(option == "both" & plen == 0)
  { stop("Opt for changepoint and order search, p.range must be provided.") }
  if(option == "cp" & plen != 0)
  { stop("Opt for changepoint search, p.range must be NULL.") }

  # set seed for reproducibility
  if(!is.null(seed)) set.seed(seed)

  if(!is.function(GA_operators$selection))  selection  = get(GA_operators$selection)
  if(!is.function(GA_operators$crossover))  crossover  = get(GA_operators$crossover)
  if(!is.function(GA_operators$mutation))   mutation   = get(GA_operators$mutation)

  ###### step 1: Initialize population
  if(class(GA_operators$population)[1] == "matrix"){
    # from input
    if(all(is.na(GA_operators$population))){stop("NA's in provided population matrix")}
    pop = GA_operators$population
  }else{
    if(!is.function(GA_operators$population)) population = get(GA_operators$population)
    # generate by function
    pop = population(popsize, p.range, N, minDist, Pchangepoint, mmax, lmax)
  }

  ## evaluate the fitness (Parallel or NOT)
  if(parallel){
    nAvaCore = detectCores()
    if(is.null(nCore)){stop(paste0("Missing number of computing cores (", nAvaCore, " cores available)."))}
    registerDoParallel(cores = nCore)
    popFit = foreach(i=1:popsize, .combine = "c") %dopar%
      (
        do.call(ObjFunc, c(list(pop[1:(pop[1,i]+plen+2),i], plen, ...)))
      )
  }else{
    popFit = rep(NA, popsize)
    for(j in 1:popsize){
      popFit[j] = do.call(ObjFunc, c(list(pop[1:(pop[1,j]+plen+2),j], plen, ...)))
      # popFit[j] = do.call(ObjFunc, c(list(pop[1:(pop[1,j]+plen+2),j], plen, Xt)))
      # popFit[j] = do.call(ObjFunc, c(list(pop[1:(pop[1,j]+plen+2),j], plen, XMat, Xt)))
      # popFit[j] = do.call(ObjFunc, c(list(pop[1:(pop[1,j]+plen+2),j], plen, XMat, Xt, logL0)))
    }
  }

  count = 0
  bestfit = rep(NA, maxgen)
  bestchrom = matrix(0, nrow=lmax, ncol=maxgen)
  repeat{
    # indicator for c("crossover", "mutation")
    #     flag[1]=1 indicating no cross-over
    #     flag[2]=1 indicating no mutation
    flag = rep(0, 2)
    ##### step 2: parents selection
    parents = selection(pop, popFit)
    dad = parents$dad
    mom = parents$mom

    ##### step 3: crossover
    a1 = runif(1)
    if(a1 <= Pcrossover){
      child = crossover(mom, dad, p.range, minDist, lmax, N)
    }else{
      child = dad
      flag[1] = 1
    }

    ##### step 4: mutation
    a2 = runif(1)
    if(a2 <= Pmutation){
      child = mutation(child, p.range, minDist, Pchangepoint, lmax, mmax, N)
    }else{
      flag[2] = 1
    }

    ##### step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum = flag[1] + flag[2]

    if (flagsum<2){
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      fitChild = do.call(ObjFunc, c(list(child[1:(child[1]+plen+2)], plen, ...)))
      # fitChild = do.call(ObjFunc, c(list(child[1:(child[1]+plen+2)], plen, XMat, Xt)))
      leastfit = max(popFit) # with largest fitness value

      if (fitChild < leastfit) {
        # indicating child is better than the worst one and replace
        pp = which.max(popFit)
        pop[,pp] = child
        popFit[pp] = fitChild
      }
    }

    # best popFit with minimized fitness value
    count = count + 1
    genbest = which.min(popFit)
    bestfit[count] = popFit[genbest]
    bestchrom[,count] = pop[,genbest]

    # check convergence once count >= maxconv
    # if after maxconv consecutive generations, the overall best does not change, then stop
    if(count >= maxconv){
      tmpbestfit = bestfit[(count-maxconv+1):count]
      decision = checkConv(tmpbestfit, maxconv, tol)
      if(monitoring){
        cat("\n My decision:", decision)
      }
      if (decision==1){
        overbestfit = bestfit[count]
        overbestchrom = bestchrom[,count]
        overbestchrom = overbestchrom[1:(overbestchrom[1]+plen+2)]
        if(monitoring){
          cat("\n==== No.", count, "Generation ====")
          cat("\n overall bestfit =", overbestfit)
          cat("\n overall bestchrom =", overbestchrom, "\n")
        }
        pp = which(is.na(bestfit))[1] - 1
        bestfit = bestfit[1:pp]
        bestchrom = bestchrom[,1:pp]
        convg = 0

        break
      }
    }

    overbestfit = bestfit[count]
    overbestchrom = bestchrom[,count]
    overbestchrom = overbestchrom[1:(overbestchrom[1]+plen+2)]
    if(monitoring){
      cat("\n==== No.", count, "Generation ====")
      cat("\n overall bestfit =", overbestfit)
      cat("\n overall bestchrom =", overbestchrom, "\n")
    }

    if (count >= maxgen){
      convg = 1
      break}
  }

  return(list(overbestfit=overbestfit, overbestchrom=overbestchrom,
              bestfit=bestfit, bestchrom=bestchrom,
              count=count, convg=convg))
}
