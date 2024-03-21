#' Genetic algorithm
#'
#' Perform the modified genetic algorithm for multiple changepoint detection.
#'
#' @param ObjFunc The fitness function to be maximized. Users can specify any R
#' functions as the fitness function with setting input as potential solution to
#' the optimization problem and returning a numerical value as the output/fitness.
#' @param n
#' @param GA_param
#' @param ga_operators
#' @param ... additional arguments that will be passed to the fitness function.
#' @Value Returns a list that has components:
#' \item{overbestfit}{Matrix of simulated presence-absence data}
#' \item{overbestchrom}{Matrix of simulated relative-abundance data}
#' \item{bestfit}{Matrix of simulated count data}
#' \item{bestchrom}{Matrix of simulated count data}
#' \item{count}{}
#' \item{convg}{}
#' @author Mo Li
#'
#' @examples
#'
#' data("throat.otu.tab")
#' otu.tab = throat.otu.tab[,colSums(throat.otu.tab>0)>1]
#'
#' fitted = MIDASim.setup(otu.tab)
#' fitted.modified = MIDASim.modify(fitted)
#' sim = MIDASim(fitted.modified, only.rel = FALSE)
#'
#' @importFrom stats runif
#' @export
GA = function(ObjFunc, n, GA_param, ga_operators, ... ){

  call = match.call()

  popsize      = GA_param$popsize
  Pcrossover   = GA_param$Pcrossover
  Pmutation    = GA_param$Pmutation
  Pchangepoint = GA_param$Pchangepoint
  minDist      = GA_param$minDist
  mmax         = GA_param$mmax
  lmax         = GA_param$lmax
  maxgen       = GA_param$maxgen
  maxconv      = GA_param$maxconv
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
  if(minDist >= Ts | minDist <= 1)
  { stop("Minimum number of locations between two changepoints invalid.") }
  if(lmax < mmax + 2 )
  { stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).") }

  # set seed for reproducibility
  if(!is.null(seed)) set.seed(seed)

  if(!is.function(ga_operators$selection))  selection  = get(ga_operators$selection)
  if(!is.function(ga_operators$crossover))  crossover  = get(ga_operators$crossover)
  if(!is.function(ga_operators$mutation))   mutation   = get(ga_operators$mutation)

  ###### step 1: Initialize population
  if(class(ga_operators$population)[1] == "matrix"){
    # from input
    if(any(is.na(population))){stop("NA's in provided population matrix")}
    population = ga_operators$population
  }else{
    if(!is.function(ga_operators$population)) population = get(ga_operators$population)
    # generate by function
    pop = population(popsize, n, minDist, Pchangepoint, mmax, lmax)
  }

  ## evaluate the fitness (Parallel or NOT)
  if(parallel){
    nAvaCore = detectCores()
    if(is.null(nCore)){stop(paste0("Missing number of computing cores (", nAvaCore, " cores available)."))}
    registerDoMC(cores = nCore)
    popFit = foreach(i=1:popsize, .combine = "c") %dopar%
      (
        do.call(ObjFunc, c(list(pop[2:(2+pop[1,i]),i], pop[1,i], ...)) )
      )
  }else{
    popFit = rep(NA, popsize)
    for(j in 1:popsize){
      m = pop[1,j]
      tau = pop[2:(2+m),j]
      popFit[j] = do.call(ObjFunc, c(list(tau, m, ...)) )
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
      child = crossover(mom, dad, minDist, lmax, n)
    }else{
      child = dad
      flag[1] = 1
    }

    ##### step 4: mutation
    a2 = runif(1)
    if(a2 <= Pmutation){
      child = mutation(minDist, Pchangepoint, lmax, mmax, n)
    }else{
      flag[2] = 1
    }

    ##### step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum = flag[1] + flag[2]

    if (flagsum<2){
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      mChild = child[1]
      tauChild = child[2:(2+mChild)]

      fitChild = do.call(ObjFunc, c(list(tauChild, mChild, ...)) )
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
      if (decision==1){
        overbestfit = bestfit[count]
        overbestchrom = bestchrom[,count]
        overbestchrom = overbestchrom[1:(overbestchrom[1]+2)]
        if(monitoring){
          cat("\n==== No.", count, "Generation ====")
          cat("\n overall bestfit =", overbestfit)
          cat("\n overall bestchrom =", overbestchrom, "\n")
        }
        pp = which(is.na(overbest))[1] - 1
        bestfit = bestfit[1:pp]
        bestchrom = bestchrom[,1:pp]
        convg = 0

        break
      }
    }

    overbestfit = bestfit[count]
    overbestchrom = bestchrom[,count]
    overbestchrom = overbestchrom[1:(overbestchrom[1]+2)]
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
