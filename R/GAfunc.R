
#--------------------------
SelectTau <- function(n, minDist, Pb, mmax){
  # this subroutine is used to select the changepoint locations with changepoint probability of Pb
  # some inputs ++++++++++++++++++
  #   n= sample size
  #   minDist= minimum locations between two changepoints
  #   Pb= the probability that a changepoint has occurred
  #   mmax= max # of changepoint
  # output ++++++++++++++++++
  #   tau= (tau_1, \ldots, tau_m)
  #   m= number of selected changepoint

  m <- 0
  tau <- 0
  i <- 1 + minDist

  repeat{
    a <- runif(1)

    if(a <= Pb){
      m = m + 1
      tau[m] = i
      i = i + minDist
      if(i >= n - minDist) {
        break
      }
    }else{
      i <- i+1
      if(i >= n - minDist){
        break
      }
    }
  }

  return(list(tau=tau, m=m))
}
# SelectTau(n=100, minDist=3, Pb=0.06, mmax=6)

#--------------------------
Initial_pop = function(popsize, n, minDist, Pb, mmax, lmax){
  # This function is used to generate initial populations for the GA
  # some inputs ++++++++++++++++++
  #   popsize = number of individual in each island
  #   Pb = the probability that a changepoint has occurred
  #   mp = minimum length of circle between two chpts
  #   mmax = max(m)
  #   m = max # of chpts
  #   lmax = max length of each chromosome
  # output ++++++++++++++++++
  #   initial population = pop

  pop = matrix(0, nrow=lmax, ncol=popsize)
  for(j in 1:popsize){
    resTau = SelectTau(n, minDist, Pb, mmax)
    tmpm = resTau$m
    tmptau = rep(NA, tmpm+1)
    if(tmpm == 0){
      tmptau[1] = n+1
    }else{
      tmptau[1:tmpm] = resTau$tau
      tmptau[tmpm+1] = n+1
    }
    pop[1,j] = tmpm            # number of changepoints
    pop[2:(2+tmpm),j] = tmptau # location of changepoints
  }

  return(pop)
}
# Initial_pop(popsize=30, n=100, minDist=3, Pb=0.06, mmax=100/2 - 1, lmax=2+(100/2 - 1))

#--------------------------
SelectTauChild = function(mom, dad, minDist, lmax, child_mmax, n){
  # This function is used to produce child
  # some inputs ++++++++++++++++++
  #   mom= selected chromosome representation with lower fitness function value
  #   dad= selected chromosome representation with larger fitness function value
  #   minDist = minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   child_mmax= maximum number of available changepoints for child
  #   n= sample size
  # outputs ++++++++++++++++++
  #   child  = the chromosome representation produced for offspring

  # 1). combine dad's and mom's changepoints
  mmom = mom[1]
  mdad = dad[1]
  parents = c(mom[2:(2+mmom-1)], dad[2:(2+mdad-1)])

  # 2). remove same/0/n+1 changepoints and sort
  parents = sort(unique(parents))
  parents = parents[!parents %in% c(0,n+1)]
  parents_count = length(parents)

  # 3). select subset from parents
  tauChild = 0
  mChild = 0
  tmpm = 1
  tmptau = parents[tmpm]

  repeat{
    a = runif(1)
    if(a > 0.5){
      mChild = mChild + 1
      tauChild[mChild] = tmptau
      tmptau = tmptau + 2*minDist
      if(tmptau >= n-minDist){break} # boundary changepoints limits
    }
    if(tmpm >= parents_count){break} # if the number of changepoints of child is larger than parents, break

    temp = parents[tmpm+1]
    if(tmptau < temp){
      tmptau = temp
    }else{
      tmpm = tmpm + 1
      if(tmpm >= parents_count){break}
      tmptau = parents[tmpm+1]
    }

    tmpm = tmpm + 1
    if(tmptau >= n-minDist){break}
  }

  tauChild[mChild+1] = n+1

  return(list(tauChild=tauChild, mChild=mChild))
}
# SelectTauChild(mom, dad, minDist, lmax, child_mmax, n)

#--------------------------
offspring = function(mom, dad, minDist, lmax, n){
  # This function is used to produce offspring
  # some inputs ++++++++++++++++++
  #   mom= selected chromosome representation with lower fitness function value
  #   dad= selected chromosome representation with larger fitness function value
  #   minDist = minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   n= sample size
  # outputs ++++++++++++++++++
  #   child  = the chromosome representation produced for offspring

  child = rep(0, lmax)

  # child number of changepoints
  child_mmax = dad[1] + mom[1]
  tauchild = rep(0, child_mmax)
  if(child_mmax == 0){
    # no changepoint for both mom and dad
    child[1] = 0
    child[2] = dad[2]
  }else{
    restauchild = SelectTauChild(mom, dad, minDist, lmax, child_mmax, n)
    mChild = restauchild$mChild
    tauChild = restauchild$tauChild
    child[1] = mChild
    child[2:(2+mChild)] = tauChild
  }

  return(child)
}
# offspring(mom, dad, minDist, lmax, n)

# #--------------------------
# mutation = function(minDist, Pb, lmax, mmax, n){
#   # This function is used to generate a new individual as the mutated child (could be customized if other parameter involved)
#   # some inputs ++++++++++++++++++
#   #   minDist= minimum distances between two adjacent changepoints
#   #   Pb= prob of changepoints for every time series
#   #   lmax= max length of chromosome
#   #   mmax= max number of changepoints
#   #   n= sample size
#   # outputs ++++++++++++++++++
#   #   childMut  = the chromosome representation produced from mutation
#
#   childMut = rep(0, lmax)
#
#   resTau = SelectTau(n, minDist, Pb, mmax)
#   mChild = resTau$m
#   tauChild = resTau$tau
#
#   tauNew = rep(0, mChild+1)
#   if(mChild==0){
#     tauNew[1] = n+1
#   }else{
#     tauNew[1:mChild] = tauChild[1:mChild]
#     tauNew[mChild+1] = n+1
#   }
#   childMut[1] = mChild
#   childMut[2:(2+mChild)] = tauNew
#
#   return(childMut)
# }
# # mutation(minDist, Pb, lmax, mmax, n)

#--------------------------
Newpopulation = function(ObjFunc1, pop, fit, popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, n, monitoring, ...){
  # This function is used to form new population
  # some inputs ++++++++++++++++++
  #   pop= population
  #   fit= fitness evaluated for population
  #   popsize= size of each island/subpop
  #   minDist= minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   mmax= max number of changepoints
  #   Pc= prob of crossover
  #   Pm= prob of mutation
  #   Pb= prob of changepoints for every time series
  #   maxgen= for each subpopulation, after maxgen then apply migration
  #   n= sample size
  #   X_hour= categorical time series
  #   XMat= Design matrix including covariate other than changepoint
  #   penalty= selection criterion to choose
  # outputs ++++++++++++++++++
  #   pop= updated population
  #   fit= updated population fitness
  #   bestfit  = currnt minimum of fitness values
  #   bestchrom = chromosome representation of the individual associated with bestfit

  count = 1
  repeat{
    # indicator for c("crossover", "mutation")
    #     flag[1]=1 indicating no cross-over
    #     flag[2]=1 indicating no mutation
    flag = rep(0, 2)

    ## step 3: select a pair of parents with linear ranking
    # the least fit (largest value with selectrank = 0) has the lowest probability
    # the best fit (smallest value with selectrank = popsize-1) has the largest probability
    selectrank = popsize - rank(fit)
    probs0 = 2.0*selectrank/(popsize*(popsize-1))
    parentI = sample(1:popsize, 2, prob=probs0)
    tmp = pop[,parentI]
    tmpp = order(fit[parentI], decreasing = FALSE)
    # dad has (better fit/smaller value/larger rank) than mom
    dad = tmp[,tmpp[1]]
    mom = tmp[,tmpp[2]]

    ## step 4-1: single point crossover
    a1 = runif(1)
    if(a1 <= Pc){
      child = offspring(mom, dad, minDist, lmax, n)
    }else{
      child = dad
      flag[1] = 1
    }

    ## step 4-2: mutation
    a2 = runif(1)
    if(a2 <= Pm){
      child = mutation(minDist, Pb, lmax, mmax, n)
    }else{
      flag[2] = 1
    }

    ## step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum = flag[1] + flag[2]
    if (flagsum<2){
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      mChild = child[1]
      tauChild = child[2:(2+mChild)]

      fitChild = ObjFunc1(tauChild, mChild)

      leastfit = max(fit) # with largest fitness value

      if (fitChild < leastfit) {
        # indicating child is better than the worst one and replace
        pp = which.max(fit)
        pop[,pp] = child
        fit[pp] = fitChild
      }
      count = count + 1
    }

    # best fit with minimized fitness value
    bestfit = min(fit)
    bestchrom = pop[, which.min(fit)]

    if(monitoring){
      cat("\n bestfit =", bestfit)
      cat("\n bestchrom =", bestchrom, "\n")
    }

    # check: after every maxgen generations, apply migration in GA.Main()
    if (count >= maxgen){break}
  }

  return(list(pop=pop, fit=fit, bestfit=bestfit, bestchrom=bestchrom))
}
# Newpopulation(pop, fit, popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, n, X_hour, XMat, penalty)

#-------------------------- Check Convergence
checkConv = function(a, maxconv, tol){
  # This function is used to check convergence criterion for overall GA
  # some inputs ++++++++++++++++++
  #   a= the best fitness values for maxconv consecutive migrations
  #   maxconv= if maxconv consecutive migrations, the overall best does not change, then stop
  #   tol= tolerance level for iterations
  # outputs ++++++++++++++++++
  #   decision= 1 means continue and 0 means stop
  i<-1
  repeat{
    diff = abs(a[i+1]-a[i])
    if (diff < tol){
      i = i+1
      if (i >= maxconv) {
        decision = 1
        break
      }
    }
    else{
      decision = 0
      break
    }
  }
  return(decision)
}
# checkConv(a, maxconv, tol)

#-------------------------- Genetic Algorithm Main Function is to minimize
GA.Main <- function(ObjFunc, n, ..., IslandGA_param){
  # GA: island models
  #     chromosome representation: real number: a vector of m, tau=(tau_1, ... tau_m, tau_{m+1}=n+1)
  #     selection: linear ranking: least fit has rank 0 and best fit has rank popsize-1
  #     crossover: apply crossover to the chosen parents with prob pc.
  #                if r1 <= pc then apply crossover and produce only one child; otherwise, child=best fif of parents
  #                p: either from mom or dad with equal probs
  #                tau: selected from all possible changepoints for mom and dad
  #     mutation: apply mutation to the child (the result of crossover) (#with rate pm=1/n)
  #               p: if r1 < 0.5, p=child's p; otherwise, select p from 0,1,2
  #               tau: if r1 < 0.5, tau=child's tau; otherwise, select new tau
  #     new generation: steady-state (replace least fit in the current pop with child if child's fit is better)
  #     migration: apply migration after maxgen generations
  #                the least fit in each subpop is replaced with one selected at random from the best fits of all subpops
  #     stopping: converge: if the overall best fit after each migration does not change for maxconv consecutive migrations
  #               ternimate: if the total number of migrations is larger than maxMig
  # some inputs ++++++++++++++++++
  #   ObjFunc= the objective function to minimize
  #   n= sample size
  #   popsize= size of each island/subpop
  #   islandsize= number of the island/subpop
  #   Pc= prob of crossover
  #   Pm= prob of mutation
  #   Pb= prob of changepoints for every time series
  #   minDist= minimum distance between two adjacent changepoints
  #   mmax= max number of changepoints
  #   lmax= max length of chromosome
  #   maxMig= max # for the migration, after maxMig migration then stop
  #   maxgen= for each subpopulation, after maxgen then apply migration
  #   maxconv= if maxconv consecutive migrations, the overall best does not change, then stop
  #   monitoring= if we want to have bestfit and best chrome displayed during GA
  #   tol = tolerance for changes in bestfits
  #   penalty= selection criterion to choose
  #   ...= additional arguments on to a subsequent function.
  # outputs ++++++++++++++++++
  #   bestfit  = min value of MDL
  #   bestchrom = optimal m, p, and tau

  cat("\n\n This island-GA is to minimize ... \n\n")

  ObjFunc1 = function(tau, m) ObjFunc(tau, m, ...)

  popsize    = IslandGA_param$popsize
  islandSize = IslandGA_param$islandSize
  Pc         = IslandGA_param$Pc
  Pm         = IslandGA_param$Pm
  Pb         = IslandGA_param$Pb
  minDist    = IslandGA_param$minDist
  mmax       = IslandGA_param$mmax
  lmax       = IslandGA_param$lmax
  maxMig     = IslandGA_param$maxMig
  maxgen     = IslandGA_param$maxgen
  maxconv    = IslandGA_param$maxconv
  monitoring = IslandGA_param$monitoring
  tol        = IslandGA_param$tol

  island = array(0, dim=c(lmax, popsize, islandSize))
  islandFit = matrix(0, nrow=popsize, ncol=islandSize)
  Bfit = rep(0, islandSize)
  Bchrom = matrix(0, nrow=lmax, ncol=islandSize)

  # step 1: generate initial population for each island
  # step 2: evaluate the fitness of the initial individuals in the population
  for(k in 1:islandSize){
    pop = Initial_pop(popsize, n, minDist, Pb, mmax, lmax)
    island[,,k] = pop
    for(j in 1:popsize){
      m = pop[1,j]
      tau = pop[2:(2+m),j]
      islandFit[j,k] = ObjFunc1(tau, m)
    }
  }

  countMig = 0
  overbest = rep(0, maxMig)
  overbestChrom = matrix(0, nrow=lmax, ncol=maxMig)
  repeat{
    # step 2,3,4,5: select parents, crossover, mutation, new pop
    for(k in 1:islandSize){
      resNewpop = Newpopulation(ObjFunc1=ObjFunc1, pop=island[,,k], fit=islandFit[,k], popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, n, monitoring=F)
      # update bestfit in each island
      Bfit[k] = resNewpop$bestfit
      Bchrom[,k] = resNewpop$bestchrom
      # update island chromosomes
      island[,,k] = resNewpop$pop
      islandFit[,k] = resNewpop$fit
    }

    # step 6: migration
    for (k in 1:islandSize){
      pleast = which.max(islandFit[,k])
      leastfit = islandFit[pleast,k]
      # replace the worst in kth island with the best from another randomly selected island
      pisland = sample(1:islandSize, 1)
      if(pisland != k){
        island[,pleast,k] = Bchrom[,pisland]
        islandFit[pleast,k] = Bfit[pisland]
      }
    }
    # update the overall bestfit
    for (k in 1:islandSize){
      pbest = which.min(islandFit[,k])
      Bchrom[,k] = island[,pbest,k]
      Bfit[k] = islandFit[pbest,k]
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
# GA.Main(ObjFunc, n, popsize, islandSize, Pc, Pm, Pb, minDist, mmax, lmax, maxMig, maxgen, maxconv, monitoring, tol, X_hour, XMat, penalty)


#-------------------------- Genetic Algorithm Main Function is to minimize
GA.cls <- function(ObjFunc, n, ..., clsGA_param){
  # GA: Classical genetic algorithm
  #     chromosome representation: real number: a vector of m, tau=(tau_1, ... tau_m, tau_{m+1}=n+1)
  #     selection: linear ranking: least fit has rank 0 and best fit has rank popsize-1
  #     crossover: apply crossover to the chosen parents with prob pc.
  #                if r1 <= pc then apply crossover and produce only one child; otherwise, child=best fif of parents
  #                p: either from mom or dad with equal probs
  #                tau: selected from all possible changepoints for mom and dad
  #     mutation: apply mutation to the child (the result of crossover) (#with rate pm=1/n)
  #               p: if r1 < 0.5, p=child's p; otherwise, select p from 0,1,2
  #               tau: if r1 < 0.5, tau=child's tau; otherwise, select new tau
  #     new generation: steady-state (replace least fit in the current pop with child if child's fit is better)
  #     stopping: converge: if the overall best fit after each generation does
  #               not change for maxconv consecutive generation ternimate;
  # some inputs ++++++++++++++++++
  #   ObjFunc= the objective function to minimize
  #   n= sample size
  #   popsize= size of each generation
  #   Pc= prob of crossover
  #   Pm= prob of mutation
  #   Pb= prob of changepoints for every time series
  #   minDist= minimum distance between two adjacent changepoints
  #   mmax= max number of changepoints
  #   lmax= max length of chromosome
  #   maxconv= if maxconv consecutive migrations, the overall best does not change, then stop
  #   monitoring= if we want to have bestfit and best chrome displayed during GA
  #   tol = tolerance for changes in bestfits
  #   ...= additional arguments on to a subsequent function.
  # outputs ++++++++++++++++++
  #   bestfit  = the optimized function value
  #   bestchrom = optimal m, p, and tau

  cat("\n\n This classical Genetic Algorithm is to minimize ... \n\n")

  ObjFunc1 = function(tau, m) ObjFunc(tau, m, ...)

  popsize    = clsGA_param$popsize
  Pc         = clsGA_param$Pc
  Pm         = clsGA_param$Pm
  Pb         = clsGA_param$Pb
  minDist    = clsGA_param$minDist
  mmax       = clsGA_param$mmax
  lmax       = clsGA_param$lmax
  maxgen     = clsGA_param$maxgen
  monitoring = clsGA_param$monitoring
  tol        = clsGA_param$tol

  # step 1: generate initial population
  populationFit = rep(0, popsize)
  population = Initial_pop(popsize, n, minDist, Pb, mmax, lmax)

  # step 2: evaluate the fitness of the initial individuals in the population
  for(j in 1:popsize){
    m = population[1,j]
    tau = population[2:(2+m),j]
    populationFit[j] = ObjFunc1(tau, m)
  }

  # step 2,3,4,5: select parents, crossover, mutation, new pop
  resNewpop = Newpopulation(ObjFunc1=ObjFunc1, pop=population, fit=populationFit, popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, n, monitoring)
  # update population chromosomes
  population = resNewpop$pop
  populationFit = resNewpop$fit

  # bestfit
  bestfit = resNewpop$bestfit
  bestchrom = resNewpop$bestchrom

  return(list(bestfit=bestfit, bestchrom=bestchrom, countMig=NULL))
}



selection_linearrank = function(pop, popFit){
  ## linear ranking
  # least fit (largest value with selectrank = 0) has the lowest probability
  # best fit (smallest value with selectrank = popsize-1) has the largest probability
  popsize = length(popFit)

  selectrank = popsize - rank(popFit)
  probs0 = 2.0*selectrank/(popsize*(popsize-1))
  parentI = sample(1:popsize, 2, prob=probs0)
  tmp = pop[,parentI]
  tmpp = order(popFit[parentI], decreasing = FALSE)
  # dad has (better fit/smaller value/larger rank) than mom
  dad = tmp[,tmpp[1]]
  mom = tmp[,tmpp[2]]

  return(list(dad=dad, mom=mom))
}

#--------------------------
offspring_uniformcrossover = function(mom, dad, minDist, lmax, n){
  # This function is used to produce offspring
  # some inputs ++++++++++++++++++
  #   mom= selected chromosome representation with lower fitness function value
  #   dad= selected chromosome representation with larger fitness function value
  #   minDist = minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   n= sample size
  # outputs ++++++++++++++++++
  #   child  = the chromosome representation produced for offspring

  child = rep(0, lmax)

  # child number of changepoints
  child_mmax = dad[1] + mom[1]
  tauchild = rep(0, child_mmax)
  if(child_mmax == 0){
    # no changepoint for both mom and dad
    child[1] = 0
    child[2] = dad[2]
  }else{
    restauchild = SelectTauChild(mom, dad, minDist, lmax, child_mmax, n)
    mChild = restauchild$mChild
    tauChild = restauchild$tauChild
    child[1] = mChild
    child[2:(2+mChild)] = tauChild
  }

  return(child)
}

#--------------------------
mutation = function(minDist, Pb, lmax, mmax, n){
  # This function is used to generate a new individual as the mutated child (could be customized if other parameter involved)
  # some inputs ++++++++++++++++++
  #   minDist= minimum distances between two adjacent changepoints
  #   Pb= prob of changepoints for every time series
  #   lmax= max length of chromosome
  #   mmax= max number of changepoints
  #   n= sample size
  # outputs ++++++++++++++++++
  #   childMut  = the chromosome representation produced from mutation

  childMut = selectTau_cpp(n, minDist, Pb, mmax, lmax)

  return(childMut)
}
# mutation(minDist, Pb, lmax, mmax, n)


# NewpopulationIsland = function(ObjFunc, selection, crossover, mutation, pop, fit, popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, n, monitoring, ...){
#   # This function is used to form new population
#   # some inputs ++++++++++++++++++
#   #   pop= population
#   #   fit= fitness evaluated for population
#   #   popsize= size of each island/subpop
#   #   minDist= minimum distances between two adjacent changepoints
#   #   lmax= max length of chromosome
#   #   mmax= max number of changepoints
#   #   Pc= prob of crossover
#   #   Pm= prob of mutation
#   #   Pb= prob of changepoints for every time series
#   #   maxgen= for each subpopulation, after maxgen then apply migration
#   #   n= sample size
#   #   X_hour= categorical time series
#   #   XMat= Design matrix including covariate other than changepoint
#   #   penalty= selection criterion to choose
#   # outputs ++++++++++++++++++
#   #   pop= updated population
#   #   fit= updated population fitness
#   #   bestfit  = currnt minimum of fitness values
#   #   bestchrom = chromosome representation of the individual associated with bestfit
#
#   count = 1
#   repeat{
#     # indicator for c("crossover", "mutation")
#     #     flag[1]=1 indicating no cross-over
#     #     flag[2]=1 indicating no mutation
#     flag = rep(0, 2)
#     ##### step 2: parents selection
#     parents = selection(pop, fit)
#     dad = parents$dad
#     mom = parents$mom
#
#     ##### step 3: crossover
#     a1 = runif(1)
#     if(a1 <= Pc){
#       child = offspring(mom, dad, minDist, lmax, n)
#     }else{
#       child = dad
#       flag[1] = 1
#     }
#
#     ## step 4-2: mutation
#     a2 = runif(1)
#     if(a2 <= Pm){
#       child = mutation(minDist, Pb, lmax, mmax, n)
#     }else{
#       flag[2] = 1
#     }
#
#     ## step 5: form new generation
#     #  steady state method:
#     #     replace the least fit in the current pop with child if child is better.
#     flagsum = flag[1] + flag[2]
#     if (flagsum<2){
#       # flagsum < 2 indicating new individual produced and fitness evaluation needed
#       mChild = child[1]
#       tauChild = child[2:(2+mChild)]
#       fitChild = do.call(ObjFunc, c(list(tauChild, mChild, ...)) )
#       leastfit = max(fit) # with largest fitness value
#
#       if (fitChild < leastfit) {
#         # indicating child is better than the worst one and replace
#         pp = which.max(fit)
#         pop[,pp] = child
#         fit[pp] = fitChild
#       }
#       count = count + 1
#     }
#
#     # best fit with minimized fitness value
#     bestfit = min(fit)
#     bestchrom = pop[, which.min(fit)]
#
#     if(monitoring){
#       cat("\n bestfit =", bestfit)
#       cat("\n bestchrom =", bestchrom, "\n")
#     }
#
#     # check: after every maxgen generations, apply migration in GA.Main()
#     if (count >= maxgen){break}
#   }
#
#   return(list(pop=pop, fit=fit, bestfit=bestfit, bestchrom=bestchrom))
# }

NewpopulationIsland = function(ObjFunc, selection, crossover, mutation, pop, fit, popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, n, ...){
  # This function is used to form new population
  # some inputs ++++++++++++++++++
  #   pop= population
  #   fit= fitness evaluated for population
  #   popsize= size of each island/subpop
  #   minDist= minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   mmax= max number of changepoints
  #   Pc= prob of crossover
  #   Pm= prob of mutation
  #   Pb= prob of changepoints for every time series
  #   maxgen= for each subpopulation, after maxgen then apply migration
  #   n= sample size
  #   X_hour= categorical time series
  #   XMat= Design matrix including covariate other than changepoint
  #   penalty= selection criterion to choose
  # outputs ++++++++++++++++++
  #   pop= updated population
  #   fit= updated population fitness
  #   bestfit  = currnt minimum of fitness values
  #   bestchrom = chromosome representation of the individual associated with bestfit

  count = 1
  repeat{
    # indicator for c("crossover", "mutation")
    #     flag[1]=1 indicating no cross-over
    #     flag[2]=1 indicating no mutation
    flag = rep(0, 2)
    ##### step 2: parents selection
    parents = selection(pop, fit)
    dad = parents$dad
    mom = parents$mom

    ##### step 3: crossover
    a1 = runif(1)
    if(a1 <= Pc){
      child = offspring(mom, dad, minDist, lmax, n)
    }else{
      child = dad
      flag[1] = 1
    }

    ## step 4-2: mutation
    a2 = runif(1)
    if(a2 <= Pm){
      child = mutation(minDist, Pb, lmax, mmax, n)
    }else{
      flag[2] = 1
    }

    ## step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum = flag[1] + flag[2]
    if (flagsum<2){
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      mChild = child[1]
      tauChild = child[2:(2+mChild)]
      fitChild = do.call(ObjFunc, c(list(tauChild, mChild, ...)) )
      # fitChild = BinSearch.BIC(tauChild, mChild, Xt)
      leastfit = max(fit) # with largest fitness value

      if (fitChild < leastfit) {
        # indicating child is better than the worst one and replace
        pp = which.max(fit)
        pop[,pp] = child
        fit[pp] = fitChild
      }
      count = count + 1
    }

    # check: after every maxgen generations, apply migration in GA.Main()
    if (count >= maxgen){break}
  }

  return(rbind(fit, pop))
}
