#' Island model based genetic algorithm
#'
#' Perform the modified island-based genetic algorithm (IslandGA) for multiple changepoint detection.
#' Minimization of an objective function using genetic algorithm (GA).
#' The algorithm can be run sequentially or in explicit parallelisation.
#' @param ObjFunc The fitness function to be minimized. Users can specify any R or Rcpp
#' function as the fitness function, setting the input as the potential solution to
#' the optimization problem and returning a numerical value as the output/fitness.
#' Depending on the user-specified chromosome representation, the optimization task
#' can be changepoint detection only or changepoint detection plus model order selection,
#' which can be specified via the \code{option} parameter. When
#' \code{option="both"}, the list \code{p.range} must be specified to give the range
#' of model orders.
#' @param N The sample size of the time series.
#' @param p.range Default is \code{NULL} for only changepoint detection. If
#' \code{p.range} is specified as a list object, which contains the range of
#' each model order parameters for order selection (integers). The number of
#' order parameters must be equal to the length of \code{p.range}.
#' @param popSize An integer representing the total number of individuals in each generation, which
#' equal to the number of islands multiplied by the size of each island (i.e., \code{popSize = number_of_islands × Islandsize}).
#' @param numIslands An integer representing the number of islands (sub-populations).
#' @param pcrossover The probability that the crossover operator applies on two individual chromosomes.
#' @param pmutation The probability that the mutation operator applies on one individual chromosome.
#' @param pchangepoint The probability that a changepoint has occurred. User could change this probability based on domain knowledge and the time series length.
#' @param minDist The minimum length between two adjacent changepoints. Default value equals to one.
#' @param mmax The maximum possible number of changepoints in the data set. For a time series of length 1000 and we only want to detect the changepoint (\code{option="cp"}), the default value is 499. The suggested value should be based on the length of the time series. For instance, if a time series has length N, the recommended \code{mmax} should be \code{N/2-1}. It is suggested to add the number of model hyperparameters if both changepoint detection and model order selection tasks are of-interested simultaneously (\code{option="both"}).
#' @param lmax The maximum possible length of the chromosome representation. For a time series of length 1000 and we only want to detect the changepoint (\code{option="cp"}), the default value is 501. The suggested value should be based on the length of the time series. For instance, if a time series has length N, the recommended \code{lmax} should be \code{2+N/2-1}. It is suggested to add the number of model hyperparameters if both changepoint detection and model order selection tasks are of-interested simultaneously (\code{option="both"}).
#' @param maxMig An integer indicates the maximum number of migrations. After conducting \code{maxMig} migrations, the island-based GA algorithm stops.
#' @param maxgen An integer indicates the maximum number of generations that each island (subpopulation) undergoes before migration. It also determines the frequency of migration. The migration, as described earlier, will occur after \code{maxgen} generations for each sub-population (island).
#' @param maxconv An integer value is also used for algorithm termination. If the overall best-fitted value remains unchanged after \code{maxconv} consecutive migrations, the island-based GA algorithm will terminate.
#' @param option A character string: either \code{"cp"} for changepoint detection only, or \code{"both"} for changepoint detection and model order selection.
#' @param monitoring Logical. If \code{TRUE}, prints intermediate GA progress.
#' @param parallel A logical value \code{TRUE} or \code{FALSE}. Default value is FALSE. It indicates whether use multiple threads to parallel compute the individual fittness function values..
#' @param nCore An integer. Default value is \code{NULL}. It represents the number of cores used in parallel computing. It must be specified if setting \code{parallel=TRUE}.
#' @param tol An numerical value. Default is \code{1e-05}. The tolerance level for deciding GA to stop.
#' @param seed An integer. Default value \code{NULL}. An single integer allows function produce reproducible results.
#' @param popInitialize A function. It should be designed for initializing a population. The default population initialization is random initialization with some imposed constraints. See \code{\link{random_population}} for example. The function returned object is a matrix, \code{pop}. The users can specified their own \code{population} function. It could also be a matrix object, which contain the user specified chromosome. By default, each column represents one individual chromosome. See \code{\link{random_population}} for details.
#' @param suggestions A list object. Default value is \code{NULL}. Each element only needs to include the better or more reasonable results of the multiple changepoint locations. Having better \code{suggestions} can help GA converges faster.
#' @param selection A function. This GA operator can help select \code{mom} and \code{dad} from current generation population, where \code{dad} is set to have better fit (smaller fitness function values). The default for selection uses the linear rank selection method. See \code{\link{selection_linearrank}} for example. The function returned object is a list contain the chromosomes for \code{mom} and \code{dad}.
#' @param crossover A function. This GA operator can apply crossover to the chosen parents to produce child for next generation with specified probability. The default for crossover uses the uniform crossover method. See \code{\link{uniformcrossover}} for details in the default crossover operator. The function returned object is a vector contain the chromosomes for \code{child}.
#' @param mutation A function. This GA operator can apply mutation to the produced child with the specified probability \code{pmutation}. See \code{\link{mutation}} for details in the default mutation operator. The function returned object is a vector contain \code{child} chromosome representation.
#' @param ... additional arguments that will be passed to the fitness function.
#' @return Return an object class \code{cptgaisl-class}. See \code{\link{cptgaisl-class}} for a more detailed description.
#' @import stats
#' @import Rcpp
#' @import foreach
#' @import doParallel
#' @import RcppArmadillo
#' @import parallel
#' @importFrom methods new
#' @importFrom utils str
#' @useDynLib changepointGA
#' @export
#' @examples
#' \donttest{
#' 
#' N = 1000
#' XMatT = matrix(1, nrow=N, ncol=1)
#' Xt = ts.sim(beta=0.5, XMat=XMatT, sigma=1, phi=0.5, theta=NULL,
#'             Delta=c(2, -2), CpLoc=c(250, 750), seed=1234)
#'
#' ## Multiple changepoint detection without model order selection
#' 
#' # without suggestions
#' GAISL.res = cptgaisl(ObjFunc=ARIMA.BIC, N=N, XMat=XMatT, Xt=Xt)
#' summary(GAISL.res)
#' plot(GAISL.res, data=Xt)
#' 
#' # with suggestions
#' suggestions = list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
#' GAISL.res = cptgaisl(ObjFunc=ARIMA.BIC, N=N, suggestions=suggestions, XMat=XMatT, Xt=Xt)
#' summary(GAISL.res)
#' plot(GAISL.res, data=Xt)
#' 
#' 
#' ## Multiple changepoint detection with model order selection
#'
#' p.range=list(ar=c(0,3), ma=c(0,3))
#'
#' # without suggestions
#' GAISL.res = cptgaisl(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, 
#'             option = "both", XMat=XMatT, Xt=Xt)
#' summary(GAISL.res)
#' plot(GAISL.res, data=Xt)
#' 
#' # with suggestions
#' suggestions = list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
#' GAISL.res = cptgaisl(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, 
#'             suggestions=suggestions, option = "both", XMat=XMatT, Xt=Xt)
#' summary(GAISL.res)
#' plot(GAISL.res, data=Xt)
#' }
cptgaisl = function(ObjFunc, 
                    N, 
                    p.range=NULL, 
                    popSize=200,
                    numIslands=5,
                    pcrossover=0.95,
                    pmutation=0.15,
                    pchangepoint=0.01,
                    minDist=1,
                    mmax=NULL,
                    lmax=NULL,
                    maxMig=1000,
                    maxgen=50,
                    maxconv=100,
                    option="cp",
                    monitoring=FALSE,
                    parallel=FALSE,
                    nCore=NULL,
                    tol=1e-05,
                    seed=NULL,
                    popInitialize = "random_population",
                    suggestions = NULL,
                    selection  = "selection_linearrank",
                    crossover  = "uniformcrossover",
                    mutation   = "mutation",
                    ... ){
  
  call = match.call()
  
  if(missing(N)){ stop("The sample size must be provided") }
  
  if (is.null(mmax)) {mmax <- floor(N/2 - 1)}
  if (is.null(lmax)) {lmax <- floor(2 + N/2 - 1)}
  
  if(missing(ObjFunc)){ stop("A fitness function must be provided") }
  if(!is.function(ObjFunc)){ stop("A fitness function must be provided") }
  
  i = NULL # add for global variables declare
  plen = length(p.range)
  
  if(pcrossover < 0   | pcrossover > 1){ stop("Probability of crossover must be between 0 and 1.") }
  if(pmutation < 0    | pmutation > 1){ stop("Probability of mutation must be between 0 and 1.") }
  if(pchangepoint < 0 | pchangepoint > 1){ stop("Probability of changepoint must be between 0 and 1.") }
  if(minDist >= N | minDist < 1){ stop("Minimum number of locations between two changepoints invalid.") }
  if(lmax < mmax + 2 ){ stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).") }
  if(option == "both" & plen == 0){ stop("Opt for changepoint and order search, p.range must be provided.") }
  if(option == "cp" & plen != 0){ stop("Opt for changepoint search, p.range must be NULL.") }
  
  # set seed for reproducibility
  if(!is.null(seed)){set.seed(seed)}
  
  if(!is.function(popInitialize)){popInitialize = get(popInitialize)}
  if(!is.function(selection)){selection = get(selection)}
  if(!is.function(crossover)){crossover = get(crossover)}
  if(!is.function(mutation)){mutation = get(mutation)}
  
  object = new("cptgaisl", 
               call=call, 
               N=N,
               p.range=p.range,
               popSize=popSize,
               numIslands=numIslands,
               Islandsize=double(),
               pcrossover=pcrossover,
               pmutation=pmutation,
               pchangepoint=pchangepoint,
               minDist=minDist,
               mmax=mmax,
               lmax=lmax,
               maxMig=maxMig,
               maxgen=maxgen,
               maxconv=maxconv,
               option=option,
               monitoring=monitoring,
               parallel=parallel,
               nCore=nCore,
               tol=tol,
               seed=seed,
               suggestions=suggestions,
               Island=array(),
               IslandFit=matrix(),
               overbestchrom=vector(),
               overbestfit=double(),
               bestfit=vector(),
               countMig=integer(),
               count=integer(),
               convg=integer())
  
  ###### step 1: generate initial population for each island
  if(!is.null(suggestions)){
    if(!is.list(suggestions)){ stop("Suggested changepoints must be a List.") }
    # 1. with suggestions
    if (any(sapply(suggestions, function(x) any(is.na(x))))) {stop("NA provided in suggestions list")}
    n.suggs = length(suggestions)
    suggestions.mat = matrix(0L, nrow=lmax, ncol=n.suggs)
    for (i in seq_along(suggestions)) {
      idx = suggestions[[i]]
      if (any(idx <= 1 | idx > N)){ stop(paste0("\n No. ", i, "th suggestion is invalid.")) }
      if (object@option == "cp"){
        if (is.null(idx)){
          suggestions.mat[1:2,i] = c(0, N+1)
        }else{
          idx = idx[idx > 1 & idx <= N]  # filter valid indices
          suggestions.mat[1:(length(idx)+2),i] = c(length(idx), idx, N+1)
        }
      }else if (object@option == "both"){
        hyper.param = sapply(p.range, function(x) if (length(x) > 0) sample(min(x):max(x), 1) else NA)
        if (is.null(idx)){
          suggestions.mat[1:(2+plen),i] = c(0, hyper.param, N+1)
        }else{
          idx = idx[idx > 1 & idx <= N]  # filter valid indices
          suggestions.mat[1:(length(idx)+plen+2),i] = c(length(idx), hyper.param, idx, N+1)
        }
      }
    }
    RemainpopSize = popSize - n.suggs
    if(RemainpopSize > 0){
      # 1.1 partial from suggestions
      if(!is.function(popInitialize)) popInitialize = get(popInitialize)
      pop = popInitialize(RemainpopSize, p.range, N, minDist, pchangepoint, mmax, lmax)  
      pop = cbind(pop, suggestions.mat)
    }else if (RemainpopSize == 0){
      # 1.2 completely from suggestions
      pop = suggestions.mat
    }else{
      # 1.3 some errors
      stop("Number of suggested chromosome must be less than specified population size.")
    }
  }else{
    # 2. completely generated by popInitialize
    if(!is.function(popInitialize)) popInitialize = get(popInitialize)
    pop = popInitialize(popSize, p.range, N, minDist, pchangepoint, mmax, lmax)
  }
  
  Islandsize = floor(popSize / numIslands)
  remainder = popSize %% numIslands
  # Shuffle individual indices
  shuffled_idx = sample(1:popSize)
  Island = array(0, dim=c(lmax, Islandsize, numIslands))
  for (k in 1:numIslands) {
    idx_start = (k - 1) * Islandsize + 1
    idx_end = k * Islandsize
    idx = shuffled_idx[idx_start:idx_end]
    Island[ , , k] = pop[,idx]
  }
  if (remainder > 0) {
    leftover_idx = shuffled_idx[(numIslands * Islandsize + 1):popSize]
    message("Warning: ", remainder, " unassigned. Consider adjusting number of islands.")
    # update population size
    popSize = numIslands * Islandsize
    object@popSize = popSize
    object@numIslands = numIslands
    object@Islandsize = Islandsize
  }
  
  object@Islandsize = Islandsize
  
  ## evaluate the fitness (Parallel or NOT)
  if(parallel){
    nAvaCore = detectCores()
    if(is.null(nCore)){stop(paste0("Missing number of computing cores (", nAvaCore, " cores available)."))}
    registerDoParallel(cores = nCore)
    IslandFit = foreach(k=1:numIslands, .combine = "cbind") %dopar%
      (
        apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, ...))))
        # apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, XMat, Xt))))
      )
  }else{
    IslandFit = matrix(0, nrow=Islandsize, ncol=numIslands)
    for(k in 1:numIslands){
      IslandFit[,k] = apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, ...))))
      # IslandFit[,k] = apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[1:(x[1]+plen+2)], plen, XMat, Xt))))
    }
  }
  
  object@Island = Island
  object@IslandFit = IslandFit
  
  countMig = 0
  
  Bfit = rep(0, numIslands)
  Bchrom = matrix(0, nrow=lmax, ncol=numIslands)
  
  bestfit = rep(0, maxMig)
  bestchrom = matrix(0, nrow=lmax, ncol=maxMig)
  
  object@bestfit = bestfit
  
  repeat{
    # step 2,3,4,5: select parents, crossover, mutation, new pop
    if(parallel){
      resNewpop = foreach(k=1:numIslands) %dopar% (
        # NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
        #                     crossover=crossover, mutation=mutation,
        #                     pop=Island[,,k], fit=IslandFit[,k],
        #                     minDist, lmax, mmax,
        #                     pcrossover, pmutation, pchangepoint,
        #                     maxgen, N, p.range, XMat, Xt)
        NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
                            crossover=crossover, mutation=mutation,
                            pop=Island[,,k], fit=IslandFit[,k],
                            minDist, lmax, mmax,
                            pcrossover, pmutation, pchangepoint,
                            maxgen, N, p.range, ...)
      )
      for(k in 1:numIslands){
        tmpfit = resNewpop[[k]][1,]
        tmppop = resNewpop[[k]][-1,]
        Bfit[k] = min(tmpfit)
        Bchrom[,k] = tmppop[, which.min(tmpfit)]
        Island[,,k] = tmppop
        IslandFit[,k] = tmpfit
      }
    }else{
      for(k in 1:numIslands){
        # resNewpop = NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
        #                                 crossover=crossover, mutation=mutation,
        #                                 pop=Island[,,k], fit=IslandFit[,k],
        #                                 minDist, lmax, mmax,
        #                                 Pc=pcrossover, Pm=pmutation, Pb=pchangepoint,
        #                                 maxgen, N, p.range, XMat, Xt)
        resNewpop = NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
                                        crossover=crossover, mutation=mutation,
                                        pop=Island[,,k], fit=IslandFit[,k],
                                        minDist, lmax, mmax,
                                        pcrossover, pmutation, pchangepoint,
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
    for (k in 1:numIslands){
      pleast = which.max(IslandFit[,k])
      leastfit = IslandFit[pleast,k]
      # replace the worst in kth island with the best from another randomly selected island
      pisland = sample(1:numIslands, 1)
      if(pisland != k){
        Island[,pleast,k] = Bchrom[,pisland]
        IslandFit[pleast,k] = Bfit[pisland]
      }
    }
    # update the overall bestfit
    for (k in 1:numIslands){
      pbest = which.min(IslandFit[,k])
      Bchrom[,k] = Island[,pbest,k]
      Bfit[k] = IslandFit[pbest,k]
    }
    
    countMig = countMig + 1
    genbest = which.min(Bfit)
    bestfit[countMig] = Bfit[genbest]
    bestchrom[,countMig] = Bchrom[,genbest]
    
    object@countMig = countMig
    object@count = countMig*maxgen
    object@Island = Island
    object@IslandFit = IslandFit
    object@bestfit = bestfit
    
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
        
        object@convg = 0
        object@overbestchrom = overbestchrom
        object@overbestfit = overbestfit
          
        break
      }
    }
    
    # step 8: check stopping if countMig > maxMig
    if(countMig >= maxMig){
      overbestfit = bestfit[countMig]
      overbestchrom = bestchrom[,countMig]
      overbestchrom = overbestchrom[1:(overbestchrom[1]+plen+2)]
      
      object@convg = 1
      object@overbestchrom = overbestchrom
      object@overbestfit = overbestfit
      
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
  
  return(object)
}
