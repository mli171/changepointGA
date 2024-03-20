#-------------------------- Genetic Algorithm Main Function is to minimize
# IslandGA = function(ObjFunc, n, IslandGA_param, IslandGA_operators, ... ){
#
#   call = match.call()
#
#   popsize      = IslandGA_param$popsize
#   Islandsize   = IslandGA_param$Islandsize
#   Pcrossover   = IslandGA_param$Pcrossover
#   Pmutation    = IslandGA_param$Pmutation
#   Pchangepoint = IslandGA_param$Pchangepoint
#   minDist      = IslandGA_param$minDist
#   mmax         = IslandGA_param$mmax
#   lmax         = IslandGA_param$lmax
#   maxMig       = IslandGA_param$maxMig
#   maxgen       = IslandGA_param$maxgen
#   maxconv      = IslandGA_param$maxconv
#   monitoring   = IslandGA_param$monitoring
#   parallel     = IslandGA_param$parallel
#   nCore        = IslandGA_param$nCore
#   tol          = IslandGA_param$tol
#   seed         = IslandGA_param$seed
#
#   if(missing(ObjFunc))
#   { stop("A fitness function must be provided") }
#   if(!is.function(ObjFunc))
#   { stop("A fitness function must be provided") }
#   if(Pcrossover < 0   | Pcrossover > 1)
#   { stop("Probability of crossover must be between 0 and 1.") }
#   if(Pmutation < 0    | Pmutation > 1)
#   { stop("Probability of mutation must be between 0 and 1.") }
#   if(Pchangepoint < 0 | Pchangepoint > 1)
#   { stop("Probability of changepoint must be between 0 and 1.") }
#   if(minDist >= Ts | minDist <= 1)
#   { stop("Minimum number of locations between two changepoints invalid.") }
#   if(lmax < mmax + 2 )
#   { stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).") }
#
#   # set seed for reproducibility
#   if(!is.null(seed)) set.seed(seed)
#
#   if(!is.function(IslandGA_operators$selection))  selection = get(IslandGA_operators$selection)
#   if(!is.function(IslandGA_operators$crossover))  crossover = get(IslandGA_operators$crossover)
#   if(!is.function(IslandGA_operators$mutation))   mutation  = get(IslandGA_operators$mutation)
#
#   Bfit = rep(0, Islandsize)
#   Bchrom = matrix(0, nrow=lmax, ncol=Islandsize)
#
#   ###### step 1: generate initial population for each island
#   if(class(IslandGA_operators$population)[1] == "matrix"){
#     # from input
#     if(any(is.na(population))){stop("NA's in provided population matrix")}
#     population = IslandGA_operators$population
#   }else{
#     if(!is.function(IslandGA_operators$population)) population = get(IslandGA_operators$population)
#     # generate by function
#     Island = array(0, dim=c(lmax, popsize, Islandsize))
#     for(k in 1:Islandsize){
#       Island[,,k] = population(popsize, n, minDist, Pchangepoint, mmax, lmax)
#     }
#   }
#
#   ## evaluate the fitness (Parallel or NOT)
#   if(parallel){
#     nAvaCore = detectCores()
#     if(is.null(nCore)){stop(paste0("Missing number of computing cores (", nAvaCore, " cores available)."))}
#     registerDoMC(cores = nCore)
#     IslandFit = foreach(k=1:Islandsize, .combine = "cbind") %dopar%
#       (
#         apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[2:(x[1]+2)], x[1], ...))))
#       )
#   }else{
#     IslandFit = matrix(0, nrow=popsize, ncol=Islandsize)
#     for(k in 1:Islandsize){
#       IslandFit[,k] = apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[2:(x[1]+2)], x[1], ...))))
#     }
#   }
#
#   countMig = 0
#   overbest = rep(0, maxMig)
#   overbestChrom = matrix(0, nrow=lmax, ncol=maxMig)
#   repeat{
#     # step 2,3,4,5: select parents, crossover, mutation, new pop
#     for(k in 1:Islandsize){
#       resNewpop = NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
#                                       crossover=crossover, mutation=mutation,
#                                       pop=Island[,,k], fit=IslandFit[,k],
#                                       popsize, minDist, lmax, mmax,
#                                       Pcrossover, Pmutation, Pchangepoint,
#                                       maxgen, n, monitoring=F, ...)
#       # update bestfit in each island
#       Bfit[k] = resNewpop$bestfit
#       Bchrom[,k] = resNewpop$bestchrom
#       # update island chromosomes
#       Island[,,k] = resNewpop$pop
#       IslandFit[,k] = resNewpop$fit
#     }
#
#     # step 6: migration
#     for (k in 1:Islandsize){
#       pleast = which.max(IslandFit[,k])
#       leastfit = IslandFit[pleast,k]
#       # replace the worst in kth island with the best from another randomly selected island
#       pisland = sample(1:Islandsize, 1)
#       if(pisland != k){
#         Island[,pleast,k] = Bchrom[,pisland]
#         IslandFit[pleast,k] = Bfit[pisland]
#       }
#     }
#     # update the overall bestfit
#     for (k in 1:Islandsize){
#       pbest = which.min(IslandFit[,k])
#       Bchrom[,k] = Island[,pbest,k]
#       Bfit[k] = IslandFit[pbest,k]
#     }
#
#     countMig = countMig + 1
#     genbest = which.min(Bfit)
#     overbest[countMig] = Bfit[genbest]
#     overbestChrom[,countMig] = Bchrom[,genbest]
#
#     # step 7: check convergence once countMig >= maxconv
#     if(countMig >= maxconv){
#       tmpoverbest = overbest[(countMig-maxconv+1):countMig]
#       decision = checkConv(tmpoverbest, maxconv, tol)
#       if (decision==1){
#         bestfit = overbest[countMig]
#         bestchrom = overbestChrom[,countMig]
#         break
#       }
#     }
#
#     # step 8: check stopping if countMig > maxMig
#     if(countMig >= maxMig){
#       bestfit = overbest[countMig]
#       bestchrom = overbestChrom[,countMig]
#       break
#     }
#
#     if(monitoring){
#       bestfit = overbest[countMig]
#       bestchrom = overbestChrom[,countMig]
#
#       cat("\n==== No.", countMig, "Migration ====")
#       cat("\n countMig =", countMig)
#       cat("\n bestfit =", bestfit)
#       cat("\n bestchrom =", bestchrom, "\n")
#     }
#   }
#
#   return(list(bestfit=bestfit, bestchrom=bestchrom, countMig=countMig))
# }

#-------------------------- Genetic Algorithm Main Function is to minimize
IslandGA = function(ObjFunc, n, IslandGA_param, IslandGA_operators, ... ){

  call = match.call()

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
  if(minDist >= Ts | minDist <= 1)
  { stop("Minimum number of locations between two changepoints invalid.") }
  if(lmax < mmax + 2 )
  { stop("Maximum length of chromosome needs to be larger than (maximum number of changepoints+2).") }

  # set seed for reproducibility
  if(!is.null(seed)) set.seed(seed)

  if(!is.function(IslandGA_operators$selection))  selection = get(IslandGA_operators$selection)
  if(!is.function(IslandGA_operators$crossover))  crossover = get(IslandGA_operators$crossover)
  if(!is.function(IslandGA_operators$mutation))   mutation  = get(IslandGA_operators$mutation)

  Bfit = rep(0, Islandsize)
  Bchrom = matrix(0, nrow=lmax, ncol=Islandsize)

  ###### step 1: generate initial population for each island
  if(class(IslandGA_operators$population)[1] == "matrix"){
    # from input
    if(any(is.na(population))){stop("NA's in provided population matrix")}
    population = IslandGA_operators$population
  }else{
    if(!is.function(IslandGA_operators$population)) population = get(IslandGA_operators$population)
    # generate by function
    Island = array(0, dim=c(lmax, popsize, Islandsize))
    for(k in 1:Islandsize){
      Island[,,k] = population(popsize, n, minDist, Pchangepoint, mmax, lmax)
    }
  }

  ## evaluate the fitness (Parallel or NOT)
  if(parallel){
    nAvaCore = detectCores()
    if(is.null(nCore)){stop(paste0("Missing number of computing cores (", nAvaCore, " cores available)."))}
    registerDoMC(cores = nCore)
    IslandFit = foreach(k=1:Islandsize, .combine = "cbind") %dopar%
      (
        apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[2:(x[1]+2)], x[1], ...))))
        # apply(Island[,,k], 2, function(x) BinSearch.BIC(x[2:(x[1]+2)], x[1], Xt))
      )
  }else{
    IslandFit = matrix(0, nrow=popsize, ncol=Islandsize)
    for(k in 1:Islandsize){
      IslandFit[,k] = apply(Island[,,k], 2, function(x) do.call(ObjFunc, c(list(x[2:(x[1]+2)], x[1], ...))))
      # IslandFit[,k] = apply(Island[,,k], 2, function(x) BinSearch.BIC(x[2:(x[1]+2)], x[1], Xt))
    }
  }

  countMig = 0
  overbest = rep(0, maxMig)
  overbestChrom = matrix(0, nrow=lmax, ncol=maxMig)
  repeat{
    # step 2,3,4,5: select parents, crossover, mutation, new pop
    if(parallel){
      resNewpop = foreach(k=1:Islandsize) %dopar% (
        NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
                            crossover=crossover, mutation=mutation,
                            pop=Island[,,k], fit=IslandFit[,k],
                            popsize, minDist, lmax, mmax,
                            Pcrossover, Pmutation, Pchangepoint,
                            maxgen, n, ...)
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
        resNewpop = NewpopulationIsland(ObjFunc=ObjFunc, selection=selection,
                                        crossover=crossover, mutation=mutation,
                                        pop=Island[,,k], fit=IslandFit[,k],
                                        popsize, minDist, lmax, mmax,
                                        Pcrossover, Pmutation, Pchangepoint,
                                        maxgen, n, ...)
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
    overbest[countMig] = Bfit[genbest]
    overbestChrom[,countMig] = Bchrom[,genbest]

    # step 7: check convergence once countMig >= maxconv
    if(countMig >= maxconv){
      tmpoverbest = overbest[(countMig-maxconv+1):countMig]
      decision = checkConv(tmpoverbest, maxconv, tol)
      if (decision==1){
        bestfit = overbest[countMig]
        bestchrom = overbestChrom[,countMig]
        break
      }
    }

    # step 8: check stopping if countMig > maxMig
    if(countMig >= maxMig){
      bestfit = overbest[countMig]
      bestchrom = overbestChrom[,countMig]
      break
    }

    if(monitoring){
      bestfit = overbest[countMig]
      bestchrom = overbestChrom[,countMig]

      cat("\n==== No.", countMig, "Migration ====")
      cat("\n countMig =", countMig)
      cat("\n bestfit =", bestfit)
      cat("\n bestchrom =", bestchrom, "\n")
    }
  }

  return(list(bestfit=bestfit, bestchrom=bestchrom, countMig=countMig))
}
