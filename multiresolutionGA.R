
library(foreach)
library(doMC)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/PopInitia.cpp")
source("R/GA.R")
source("R/IslandGA.R")
source("R/GAfunc.R")

cptDist = function(tau1, tau2, n){

  m = length(tau1)
  k = length(tau2)

  if(is.null(tau1) & is.null(tau2)){
    dist.res = NA
    warning("Both configurations are NULL.")
  }else{
    if(is.null(tau1) | is.null(tau2)){
      ACC = 0
    }else{
      if(any(tau1 < 0 | tau1 > n)){stop("First changepoint configuration invalid.")}
      if(any(tau2 < 0 | tau2 > n)){stop("Second changepoint configuration invalid.")}

      costs = matrix(0, nrow=m, ncol=k)
      for(i in 1:m){
        for(j in 1:k){
          costs[i,j] = abs(tau1[i] - tau2[j])/n
        }
      }
      if(m > k){costs=t(costs)}
      y = clue::solve_LSAP(costs, maximum = FALSE)
      ACC = sum(costs[cbind(seq_along(y), y)])
    }
    dist.res = abs(m-k) + ACC
  }

  return(dist.res)
}

##### Simulation for stationary AR1 process
ar1.sim = function(theta, sigma, phi, XMat, CpConfig, seed=NULL){

  m = sum(CpConfig)

  if(m==0){
    # no changepoint
    DesignX = XMat
    mst = DesignX%*%theta
  }else{
    # changepoints
    tauClc = which(CpConfig==1)
    tmptau = unique(c(tauClc, Ts+1))
    CpMat = matrix(0, nrow=Ts, ncol=length(tmptau)-1)
    for(i in 1:NCOL(CpMat)){CpMat[tmptau[i]:(tmptau[i+1]-1),i] = 1}
    DesignX = cbind(XMat, CpMat)
    mst = DesignX%*%theta
  }

  if(!is.null(seed)){set.seed(seed)}

  et = arima.sim(n=Ts, list(ar=phi), sd=sigma) # sd argument is for WN
  Z = et + mst

  attr(Z, 'DesignX') = DesignX

  return(Z)
}

thetaT = c(0.5, 1, -1, 1, -1) # intercept, cp1, cp2
sigmaT = 1
phiT = 0.75

Ts = 30*365*24
tau.prop = c(1/4, 1/2, 2.5/4, 3/4)
tauT = floor(Ts*tau.prop)
CpConfigT = rep(0, Ts)
CpConfigT[tauT] = 1

XMatT = matrix(1, nrow=Ts, ncol=1)
colnames(XMatT) = "alpha0"

Z.sim = ar1.sim(theta=thetaT, sigma=sigmaT, phi=phiT, XMat=XMatT, CpConfig=CpConfigT, seed=NULL)
# DesignXT = attributes(Z.sim)$DesignX


## objective function for changepoint search BIC
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


tim1 = Sys.time()

ga_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank_cpp",
                    crossover  = "offspring_uniformcrossover_cpp",
                    mutation   = "mutation")

# ellaClc = c(10, 8, 6, 4, 3, 2)

ellaClc = c(365*24, 30*24, 24, 10, 8, 6, 4, 3, 2)

multids.pop = matrix(NA, nrow=2+Ts/2-1, ncol=length(ellaClc))

for(ii in 1:length(ellaClc)){

  # downsampling a TS
  ella = ellaClc[ii]
  p.seq = seq(from=1, to=Ts, by=ella)
  Z.sim.ds = Z.sim[p.seq]
  Ts.ds = length(Z.sim.ds)

  GA_param = list(
    popsize      = 200,
    Pcrossover   = 0.95,
    Pmutation    = 0.15,
    Pchangepoint = 0.06,
    minDist      = 2,
    mmax         = Ts.ds/2 - 1,
    lmax         = 2 + Ts.ds/2 - 1,
    maxgen       = 10000,
    maxconv      = 5000,
    monitoring   = FALSE,
    parallel     = FALSE,
    nCore        = NULL,
    tol          = 1e-5,
    seed         = NULL
  )

  tmp1 = GA(BinSearch.BIC, n=Ts.ds, GA_param, ga_operators, Xt=Z.sim.ds)

  m.ds = tmp1$overbestchrom[1]
  if(m.ds>0){
    tau.prop.ds = tmp1$overbestchrom[2:(1+m.ds)]/Ts.ds
    tau.ds = c(m.ds, floor(Ts*tau.prop.ds), Ts+1)
  }else{
    tau.prop.ds = NULL
    tau.ds = c(0, Ts+1)
  }
  multids.pop[1:(m.ds+2),ii] = tau.ds

  cat("\n ella = ", ella, "completed!")

}

pop.mat = random_population_cpp(popsize=200, n=Ts, minDist=2, Pb=0.06, mmax=Ts/2 - 1, lmax=2 + Ts/2 - 1)

pop.mat[,1:length(ellaClc)] = multids.pop

ga_operators = list(population = pop.mat,
                    selection  = "selection_linearrank_cpp",
                    crossover  = "offspring_uniformcrossover_cpp",
                    mutation   = "mutation")
GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.06,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxgen       = 5000,
  maxconv      = 1000,
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
tmp1 = GA(BinSearch.BIC, n=Ts, GA_param, ga_operators, Xt=Z.sim)

tim2 = Sys.time()


tim3 = Sys.time()
ga_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank_cpp",
                    crossover  = "offspring_uniformcrossover_cpp",
                    mutation   = "mutation")
GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.06,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxgen       = 10000,
  maxconv      = 5000,
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
tmp2 = GA(BinSearch.BIC, n=Ts, GA_param, ga_operators, Xt=Z.sim)
tim4 = Sys.time()


dist1 = cptDist(tauT, tmp1$overbestchrom[2:(tmp1$overbestchrom[1]+1)], Ts)
dist2 = cptDist(tauT, tmp2$overbestchrom[2:(tmp2$overbestchrom[1]+1)], Ts)

tmp1$overbestfit
tmp1$overbestchrom
tmp1$count
dist1
tmp2$overbestfit
tmp2$overbestchrom
tmp2$count
dist2

tim2 - tim1
tim4 - tim3


