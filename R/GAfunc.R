
#--------------------------
SelectTau <- function(N, minDist, Pb, mmax){
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
      if(i >= N - minDist) {
        break
      }
    }else{
      i <- i+1
      if(i >= N - minDist){
        break
      }
    }
  }

  return(list(tau=tau, m=m))
}
# SelectTau(n=100, minDist=3, Pb=0.06, mmax=6)

#--------------------------
Initial_pop = function(popsize, N, minDist, Pb, mmax, lmax){
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
    resTau = SelectTau(N, minDist, Pb, mmax)
    tmpm = resTau$m
    tmptau = rep(NA, tmpm+1)
    if(tmpm == 0){
      tmptau[1] = N+1
    }else{
      tmptau[1:tmpm] = resTau$tau
      tmptau[tmpm+1] = N+1
    }
    pop[1,j] = tmpm            # number of changepoints
    pop[2:(2+tmpm),j] = tmptau # location of changepoints
  }

  return(pop)
}
# Initial_pop(popsize=30, N=100, minDist=3, Pb=0.06, mmax=100/2 - 1, lmax=2+(100/2 - 1))

#--------------------------
SelectTauChild = function(mom, dad, minDist, lmax, child_mmax, N){
  # This function is used to produce child
  # some inputs ++++++++++++++++++
  #   mom= selected chromosome representation with lower fitness function value
  #   dad= selected chromosome representation with larger fitness function value
  #   minDist = minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   child_mmax= maximum number of available changepoints for child
  #   N= sample size
  # outputs ++++++++++++++++++
  #   child  = the chromosome representation produced for offspring

  # 1). combine dad's and mom's changepoints
  mmom = mom[1]
  mdad = dad[1]
  parents = c(mom[2:(2+mmom-1)], dad[2:(2+mdad-1)])

  # 2). remove same/0/N+1 changepoints and sort
  parents = sort(unique(parents))
  parents = parents[!parents %in% c(0,N+1)]
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
      if(tmptau >= N-minDist){break} # boundary changepoints limits
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
    if(tmptau >= N-minDist){break}
  }

  tauChild[mChild+1] = N+1

  return(list(tauChild=tauChild, mChild=mChild))
}
# SelectTauChild(mom, dad, minDist, lmax, child_mmax, N)

#--------------------------
offspring = function(mom, dad, minDist, lmax, N){
  # This function is used to produce offspring
  # some inputs ++++++++++++++++++
  #   mom= selected chromosome representation with lower fitness function value
  #   dad= selected chromosome representation with larger fitness function value
  #   minDist = minimum distances between two adjacent changepoints
  #   lmax= max length of chromosome
  #   N= sample size
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
    restauchild = SelectTauChild(mom, dad, minDist, lmax, child_mmax, N)
    mChild = restauchild$mChild
    tauChild = restauchild$tauChild
    child[1] = mChild
    child[2:(2+mChild)] = tauChild
  }

  return(child)
}

#--------------------------
Newpopulation = function(ObjFunc1, pop, fit, popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, N, monitoring, ...){
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
  #   N= sample size
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
      child = offspring(mom, dad, minDist, lmax, N)
    }else{
      child = dad
      flag[1] = 1
    }

    ## step 4-2: mutation
    a2 = runif(1)
    if(a2 <= Pm){
      child = mutation(minDist, Pb, lmax, mmax, N)
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

#--------------------------
offspring_uniformcrossover = function(mom, dad, minDist, lmax, N){
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
    restauchild = SelectTauChild(mom, dad, minDist, lmax, child_mmax, N)
    mChild = restauchild$mChild
    tauChild = restauchild$tauChild
    child[1] = mChild
    child[2:(2+mChild)] = tauChild
  }

  return(child)
}

#--------------------------
mutation = function(child, p.range=NULL, minDist, Pb, lmax, mmax, N){
  # This function is used to generate a new individual as the mutated child (could be customized if other parameter involved)
  # some inputs ++++++++++++++++++
  #   prange= The range of the model order
  #   minDist= minimum distances between two adjacent changepoints
  #   Pb= prob of changepoints for every time series
  #   lmax= max length of chromosome
  #   mmax= max number of changepoints
  #   N= sample size
  # outputs ++++++++++++++++++
  #   childMut  = the chromosome representation produced from mutation

  plen = length(p.range)

  if(plen>0){
    childMut = matrix(0, nrow=lmax, 1)
    a1 = runif(1)
    if(a1 > 0.5){
      # 1. order from child
      childMut[2:(plen+1),] = child[2:(plen+1)]
      # 1.1 cpt from new
      tmpchildMut = selectTau_cpp(N=N, prange=NULL, minDist=minDist, Pb=Pb,
                                  mmax=mmax, lmax=lmax)
      childMut[1,] = tmpchildMut[1]
      childMut[(plen+2):(plen+tmpchildMut[1]+2),] = tmpchildMut[2:(tmpchildMut[1]+2)]
    }else{
      # 2. order from new
      new.p.range = rep(NA, plen)
      for(ii in 1:plen){
        tmp.p.range = setdiff(p.range[[ii]][1]:p.range[[ii]][2], child[2+ii-1])
        new.p.range[ii] = sample(tmp.p.range, 1)
      }
      childMut[2:(plen+1),] = new.p.range
      a2 = runif(1)
      if(a2 > 0.5){
        # 2.1 cpt from new
        tmpchildMut = selectTau_cpp(N=N, prange=NULL, minDist=minDist, Pb=Pb,
                                    mmax=mmax, lmax=lmax)
        childMut[1,] = tmpchildMut[1]
        childMut[(plen+2):(plen+tmpchildMut[1]+2),] = tmpchildMut[2:(tmpchildMut[1]+2)]
      }else{
        # 2.2 cpt from child
        childMut[1,] = child[1]
        childMut[(plen+2):(plen+child[1]+2),] = child[(plen+2):(plen+child[1]+2)]
      }
    }
  }else{
    tmpchildMut = selectTau_cpp(N=N, prange=NULL, minDist=minDist, Pb=Pb,
                                mmax=mmax, lmax=lmax)
    childMut = tmpchildMut
  }

  return(childMut)
}

NewpopulationIsland = function(ObjFunc, selection, crossover, mutation, pop, fit, popsize, minDist, lmax, mmax, Pc, Pm, Pb, maxgen, N, p.range, ...){
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
  #   N= sample size
  #   X_hour= categorical time series
  #   XMat= Design matrix including covariate other than changepoint
  #   penalty= selection criterion to choose
  # outputs ++++++++++++++++++
  #   pop= updated population
  #   fit= updated population fitness
  #   bestfit  = currnt minimum of fitness values
  #   bestchrom = chromosome representation of the individual associated with bestfit

  plen = length(p.range)

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
      child = crossover(mom, dad, p.range, minDist, lmax, N)
    }else{
      child = dad
      flag[1] = 1
    }

    ## step 4-2: mutation
    a2 = runif(1)
    if(a2 <= Pm){
      child = mutation(child, p.range, minDist, Pb, lmax, mmax, N)
    }else{
      flag[2] = 1
    }

    ## step 5: form new generation
    #  steady state method:
    #     replace the least fit in the current pop with child if child is better.
    flagsum = flag[1] + flag[2]
    if (flagsum<2){
      # flagsum < 2 indicating new individual produced and fitness evaluation needed
      fitChild = do.call(ObjFunc, c(list(child[1:(child[1]+plen+2)], plen, ...)))
      # fitChild = do.call(ObjFunc, c(list(child[1:(child[1]+plen+2)], plen, XMat, Xt)))
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
