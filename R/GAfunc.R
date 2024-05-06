
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

#' The example objective function for changepoint search in First order Gaussian
#' autoregressions (AR(1)) via Bayesian Information Criterion (BIC)
#'
#' This is a function to calculate
#' @param tau The provided changepoint locations (Only locations are included).
#' @param m The provided number of changepoint locations.
#' @param Xt The simulated AR(1) time series from \code{ts.sim} function.
#' @return Returned the value of the obejctive function (i.e. BIC).
#' @import stats
#' @useDynLib changepointGA
#' @export
#' @examples
#' library(changepointGA)
#'
#' Ts = 1000
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#' DeltaT = c(2, -2)
#' sigmaT = 1
#' thetaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' myts = ts.sim(theta=thetaT, XMat=XMatT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
#'
#'# candidate changepoint configuration
#' tautest = c(250, 500, 750)
#' BinSearch.BIC(tau=tautest, m=3, Xt=myts)
BinSearch.BIC = function(tau, m, Xt){

  N = length(Xt) #length of the series

  if(m == 0){
    ##Case 1, Zero Changepoint
    mu.hat = mean(Xt)
    phi.hat = sum((Xt-mu.hat)[-N]*(Xt-mu.hat)[-1])/sum((Xt-mu.hat)[-1]^2)
    Xt.hat = c(mu.hat, mu.hat + phi.hat*(Xt[-N]-mu.hat))
    sigma.hatsq = sum( (Xt-Xt.hat)^2 )/N
    BIC.obj = N*log(sigma.hatsq)+ 3*log(N) #6 always there
  }else{
    tau = tau[tau>1 & tau<N+1] #keep CPT locations only
    tau.ext = c(1,tau,(N+1)) #include CPT boundary 1 and N+1

    ## Split Xt to regimes/segments to
    ## compute phi.hat and sigma.hat.sq
    seg.len = diff(tau.ext) #length of each segments
    ff = rep(0:m, times=seg.len) ##create factors for segmentation
    Xseg = split(Xt, ff) ##Segmentation list
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    mu.hat = rep(mu.seg, seg.len)
    phi.hat = sum((Xt-mu.hat)[-N]*(Xt-mu.hat)[-1])/sum((Xt-mu.hat)[-1]^2)
    Xt.hat = c(mu.hat[1], mu.hat[-1] + phi.hat*(Xt[-N]-mu.hat[-N]))
    sigma.hatsq = sum( (Xt-Xt.hat)^2 )/N
    BIC.obj = N*log(sigma.hatsq) + (2*m + 3)*log(N)
  }
  return(BIC.obj)
}
