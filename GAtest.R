
library(foreach)
library(doMC)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/PopInitia.cpp")
source("R/GA.R")
source("R/IslandGA.R")
source("R/GAfunc.R")

##### Simulation for stationary AR1 process
ar1.sim = function(theta, sigma, phi, XMat, CpConfig, seed=NULL){

  len_par = length(par)

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

thetaT = c(0.5, 2, -2) # intercept, cp1, cp2
sigmaT = 1
phiT = 0.5

Ts = 50000
tau.prop = c(1/4, 3/4)
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
    tau = tau[tau>1 & tau<Ts+1] #keep CPT locations only
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

GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.06,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxgen       = 10000,
  maxconv      = 10000,
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)

ga_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank",
                    crossover  = "offspring_uniformcrossover",
                    mutation   = "mutation")

tim1 = Sys.time()
tmp1 = GA(BinSearch.BIC, n=Ts, GA_param, ga_operators, Xt=Z.sim)
tim2 = Sys.time()


IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.06,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 100,
  monitoring   = FALSE,
  parallel     = FALSE, ###
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)


IslandGA_operators = list(population = "random_population_cpp",
                          selection  = "selection_linearrank",
                          crossover  = "offspring_uniformcrossover",
                          mutation   = "mutation")

tim7 = Sys.time()
tmp4 = IslandGA(BinSearch.BIC, n=Ts, IslandGA_param, IslandGA_operators, Xt=Z.sim)
tim8 = Sys.time()


IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.06,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 100,
  monitoring   = FALSE,
  parallel     = TRUE, ###
  nCore        = 5,
  tol          = 1e-5,
  seed         = NULL
)


IslandGA_operators = list(population = "random_population_cpp",
                          selection  = "selection_linearrank",
                          crossover  = "offspring_uniformcrossover",
                          mutation   = "mutation")

tim9 = Sys.time()
tmp5 = IslandGA(BinSearch.BIC, n=Ts, IslandGA_param, IslandGA_operators, Xt=Z.sim)
tim10 = Sys.time()


tim2 - tim1
# tim4 - tim3
# tim6 - tim5
tim8 - tim7
tim10 - tim9

tmp1$overbestfit
# tmp2$overbestfit
# tmp3$overbestfit
tmp4$bestfit
tmp5$bestfit

# GA
tmp1$overbestchrom
# tmp2$overbestchrom
# tmp3$overbestchrom
# island GA
tmp4$bestchrom[tmp4$bestchrom!=0]
tmp5$bestchrom[tmp5$bestchrom!=0]
