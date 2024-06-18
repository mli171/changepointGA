#' Island model based genetic algorithm
#'
#' Perform the modified island-based genetic algorithm (IslandGA) for multiple changepoint detection.
#' Minimization of an objective function using genetic algorithm (GA).
#' The algorithm can be run sequentially or in explicit parallelisation.
#' @param ObjFunc The fitness function to be minimized. Users can specify any R or Rcpp
#' functions as the fitness function with setting input as potential solution to
#' the optimization problem and returning a numerical value as the output/fitness.
#'
#' @param N The sample size of the time series.
#' @param IslandGA_param A list contains the hyper-parameters for IslandGA.
#' See \code{\link{IslandGA_param}} for the details.
#' @param IslandGA_operators A list includes the functions for population initialization,
#' new individual selection, and genetic operator of crossover and mutation.
#' See \code{\link{IslandGA_operators}} for the details.
#' @param p.range Default is \code{NULL} for only changepoint detection. If
#' \code{p.range} is specified as a list object, which contains the range of
#' each model order parameters for order selection (integers). The number of
#' order parameters must be equal to the length of \code{p.range}.
#' @param ... additional arguments that will be passed to the fitness function.
#' @return Returns a list that has the following components.
#' \item{overbestfit}{The obtained minimum value of the objective function after
#' the final iteration.}
#' \item{overbestchrom}{The locations of the detected changepoints associating
#' with the \code{overbestfit} the after the final iteration.}
#' \item{bestfit}{The minimized fitness function values at each iteration.}
#' \item{bestchrom}{The detected changepoints at each iteration.}
#' \item{countMig}{The number of migrations undertaken by the IslandGA.}
#' \item{count}{The number of iterations (generations) undertaken by the island
#' genetic algorithm model.}
#' \item{convg}{An integer code indicate convergence.}
#' \itemize{
#'  \item{0} indicates the algorithm successful completion.
#'  \item{1} indicates the the total number of migrations exceeds the
#'  pre-specified \item{\code{maxMig}} limit.
#'  }
#' @import stats
#' @import Rcpp
#' @import foreach
#' @import doMC
#' @import RcppArmadillo
#' @import parallel
#' @useDynLib changepointGA
#' @export
#' @examples
#' Ts = 1000
#' betaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' sigmaT = 1
#' phiT = c(0.5, -0.5)
#' thetaT = c(0.8)
#' DeltaT = c(2, -2)
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#'
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT,
#'               Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#' TsPlotCheck(myts, tau=CpLocT)
#'
#' IslandGA_param = list(
#'   popsize      = 40,
#'   Islandsize   = 5,
#'   Pcrossover   = 0.95,
#'   Pmutation    = 0.15,
#'   Pchangepoint = 10/Ts,
#'   minDist      = 2,
#'   mmax         = Ts/2 - 1,
#'   lmax         = 2 + Ts/2 - 1,
#'   maxMig       = 500,
#'   maxgen       = 100,
#'   maxconv      = 100,
#'   option       = "cp",
#'   monitoring   = FALSE,
#'   parallel     = FALSE,
#'   nCore        = NULL,
#'   tol          = 1e-5,
#'   seed         = NULL
#' )
#' IslandGA_operators = list(population = "random_population_cpp",
#'                           selection  = "selection_linearrank_cpp",
#'                           crossover  = "offspring_uniformcrossover_cpp",
#'                           mutation   = "mutation")
#'
#' IslandGA.res = IslandGA(ObjFunc=BinSearch.BIC, N=Ts, IslandGA_param, IslandGA_operators, Xt=myts)
#' # IslandGA.res$overbestfit
#' # IslandGA.res$overbestchrom
#-------------------------- Genetic Algorithm Main Function is to minimize
IslandGA = function(ObjFunc, N, IslandGA_param, IslandGA_operators, p.range=NULL, ... ){

  call = match.call()
  plen = length(p.range)

  popsize      = IslandGA_param$popsize
  Islandsize   = IslandGA_param$Islandsize
  Pcrossover   = IslandGA_param$Pcrossover
  Pmutation    = IslandGA_param$Pmutation
  Pchangepoint = IslandGA_param$Pchangepoint
  minDist      = IslandGA_param$minDist
  mmax         = IslandGA_param$mmax
  lmax         = IslandGA_param$lmax
  maxMig       = IslandGA_param$maxMig
  maxgen       = IslandGA_param$maxgen
  maxconv      = IslandGA_param$maxconv
  option       = IslandGA_param$option
  monitoring   = IslandGA_param$monitoring
  parallel     = IslandGA_param$parallel
  nCore        = IslandGA_param$nCore
  tol          = IslandGA_param$tol
  seed         = IslandGA_param$seed

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
  if(minDist >= N | minDist <= 1)
  { stop("Minimum number of locations between two changepoints invalid.") }
  if(lmax < mmax + 2 )
  { stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).") }
  if(option == "both" & plen == 0)
  { stop("Opt for changepoint and order search, p.range must be provided.") }
  if(option == "cp" & plen != 0)
  { stop("Opt for changepoint search, p.range must be NULL.") }

  # set seed for reproducibility
  if(!is.null(seed)) set.seed(seed)

  if(!is.function(IslandGA_operators$selection))  selection = get(IslandGA_operators$selection)
  if(!is.function(IslandGA_operators$crossover))  crossover = get(IslandGA_operators$crossover)
  if(!is.function(IslandGA_operators$mutation))   mutation  = get(IslandGA_operators$mutation)

  Bfit = rep(0, Islandsize)
  Bchrom = matrix(0, nrow=lmax, ncol=Islandsize)

  ###### step 1: generate initial population for each island
  if(class(IslandGA_operators$population)[1] == "array"){
    # from input
    if(all(is.na(IslandGA_operators$population))){stop("NA's in provided population array")}
    Island = IslandGA_operators$population
  }else{
    if(!is.function(IslandGA_operators$population)) population = get(IslandGA_operators$population)
    # generate by function
    Island = array(0, dim=c(lmax, popsize, Islandsize))
    for(k in 1:Islandsize){
      Island[,,k] = population(popsize, p.range, N, minDist, Pchangepoint, mmax, lmax)
    }
  }

  ## evaluate the fitness (Parallel or NOT)
  if(parallel){
    nAvaCore = detectCores()
    if(is.null(nCore)){stop(paste0("Missing number of computing cores (", nAvaCore, " cores available)."))}
    registerDoMC(cores = nCore)
    IslandFit = foreach(k=1:Islandsize, .combine = "cbind") %dopar%
      (
        apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, ...))))
        # apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, XMat, Xt))))
      )
  }else{
    IslandFit = matrix(0, nrow=popsize, ncol=Islandsize)
    for(k in 1:Islandsize){
      IslandFit[,k] = apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, ...))))
      # IslandFit[,k] = apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, XMat, Xt))))
    }
  }

  countMig = 0
  bestfit = rep(0, maxMig)
  bestchrom = matrix(0, nrow=lmax, ncol=maxMig)
  repeat{
    # step 2,3,4,5: select parents, crossover, mutation, new pop
    if(parallel){
      resNewpop = foreach(k=1:Islandsize) %dopar% (
        # NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
        #                     crossover=crossover, mutation=mutation,
        #                     pop=Island[,,k], fit=IslandFit[,k],
        #                     popsize, minDist, lmax, mmax,
        #                     Pcrossover, Pmutation, Pchangepoint,
        #                     maxgen, N, p.range, XMat, Xt)
        NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
                            crossover=crossover, mutation=mutation,
                            pop=Island[,,k], fit=IslandFit[,k],
                            popsize, minDist, lmax, mmax,
                            Pcrossover, Pmutation, Pchangepoint,
                            maxgen, N, p.range, ...)
      )
      for(k in 1:Islandsize){
        tmpfit = resNewpop[[k]][1,]
        tmppop = resNewpop[[k]][-1,]
        Bfit[k] = min(tmpfit)
        Bchrom[,k] = tmppop[, which.min(tmpfit)]
        Island[,,k] = tmppop
        IslandFit[,k] = tmpfit
      }
    }else{
      for(k in 1:Islandsize){
        # resNewpop = NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
        #                                 crossover=crossover, mutation=mutation,
        #                                 pop=Island[,,k], fit=IslandFit[,k],
        #                                 popsize, minDist, lmax, mmax,
        #                                 Pc=Pcrossover, Pm=Pmutation, Pb=Pchangepoint,
        #                                 maxgen, N, p.range, XMat, Xt)
        resNewpop = NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
                                        crossover=crossover, mutation=mutation,
                                        pop=Island[,,k], fit=IslandFit[,k],
                                        popsize, minDist, lmax, mmax,
                                        Pcrossover, Pmutation, Pchangepoint,
                                        maxgen, N, p.range, ...)
        tmpfit = resNewpop[1,]
        tmppop = resNewpop[-1,]
        # update bestfit in each island
        Bfit[k] = min(tmpfit)
        Bchrom[,k] = tmppop[, which.min(tmpfit)]
        # update island chromosomes
        Island[,,k] = tmppop
        IslandFit[,k] = tmpfit
      }
    }

    # step 6: migration
    for (k in 1:Islandsize){
      pleast = which.max(IslandFit[,k])
      leastfit = IslandFit[pleast,k]
      # replace the worst in kth island with the best from another randomly selected island
      pisland = sample(1:Islandsize, 1)
      if(pisland != k){
        Island[,pleast,k] = Bchrom[,pisland]
        IslandFit[pleast,k] = Bfit[pisland]
      }
    }
    # update the overall bestfit
    for (k in 1:Islandsize){
      pbest = which.min(IslandFit[,k])
      Bchrom[,k] = Island[,pbest,k]
      Bfit[k] = IslandFit[pbest,k]
    }

    countMig = countMig + 1
    genbest = which.min(Bfit)
    bestfit[countMig] = Bfit[genbest]
    bestchrom[,countMig] = Bchrom[,genbest]

    # step 7: check convergence once countMig >= maxconv
    if(countMig >= maxconv){
      tmpbest = bestfit[(countMig-maxconv+1):countMig]
      decision = checkConv(tmpbest, maxconv, tol)
      if(monitoring){
        cat("\n My decision:", decision)
      }
      if (decision==1){
        overbestfit = bestfit[countMig]
        overbestchrom = bestchrom[,countMig]
        overbestchrom = overbestchrom[1:(overbestchrom[1]+plen+2)]
        convg = 0

        break
      }
    }

    # step 8: check stopping if countMig > maxMig
    if(countMig >= maxMig){
      overbestfit = bestfit[countMig]
      overbestchrom = bestchrom[,countMig]
      overbestchrom = overbestchrom[1:(overbestchrom[1]+plen+2)]
      convg = 1

      break
    }

    if(monitoring){
      overbestfit = bestfit[countMig]
      overbestchrom = bestchrom[,countMig]
      overbestchrom = overbestchrom[1:(overbestchrom[1]+plen+2)]
      cat("\n==== No.", countMig, "Migration ====")
      cat("\n countMig =", countMig)
      cat("\n overbestfit =", overbestfit)
      cat("\n overbestchrom =", overbestchrom, "\n")
    }
  }

  return(list(overbestfit=overbestfit, overbestchrom=overbestchrom,
              bestfit=bestfit, bestchrom=bestchrom,
              countMig=countMig, count=countMig*maxgen, convg=convg))
}
