
library(foreach)
library(doMC)

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

Ts = 200
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


######################################################
#------------------- GA Setting ---------------------#
######################################################

popsize    = 40         # number of individuals in each island
islandSize = 5          # number of island
Pc         = 0.95       # prob of cross-over
Pm         = 0.15       # prob of mutation
Pb         = 0.06       # prob of changepoints for every time series
minDist    = 2          # minimum number of locations between two changepoints
mmax       = Ts/2 - 1   # maximum number of changepoints in each chromosome
lmax       = 2+mmax     # max length of each chromosome
maxMig     = 1000       # maximum migration times
maxgen     = 20         # for each subpopulation, after maxgen then apply migration
maxconv    = 20         # number of consecutive migrations
monitoring = TRUE       # whether print out the output from each Migration from GA
tol        = 1e-5       # tolerance level for iterations

IslandGA_param = list(
  popsize    = popsize,
  islandSize = islandSize,
  Pc         = Pc,
  Pm         = Pm,
  Pb         = Pb,
  minDist    = minDist,
  mmax       = mmax,
  lmax       = lmax,
  maxMig     = maxMig,
  maxgen     = maxgen,
  maxconv    = maxconv,
  monitoring = F,
  tol        = tol
)


source("R/GAfunc.R")
tim1 = Sys.time()
tmp1 = GA.Main(ObjFunc=BinSearch.BIC, n=Ts, IslandGA_param=IslandGA_param, Xt=Z.sim)
tim2 = Sys.time()



GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.06,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxgen       = 1000,
  maxconv      = 1000,
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)

ga_operators = list(population = "Initial_pop",
                    selection  = "selection",
                    crossover  = "offspring",
                    mutation   = "mutation")
sourceCpp("src/PopInitia.cpp")
ga_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank",
                    crossover  = "offspring_uniformcrossover",
                    mutation   = "mutation_new")

tim3 = Sys.time()
tmp2 = GA(BinSearch.BIC, n=Ts, GA_param, ga_operators, Xt=Z.sim)
tim4 = Sys.time()



##BIC
BinSearch.BIC = function(loc.ind, Xt=xt){
  loc.ind[1]=0
  N = length(Xt) #length of the series
  m = sum(loc.ind) #Number of CPTs

  if(m == 0){
    ##Case 1, Zero Changepoint
    mu.hat = mean(Xt)
    phi.hat = sum((Xt-mu.hat)[-N]*(Xt-mu.hat)[-1])/sum((Xt-mu.hat)[-1]^2)
    Xt.hat = c(mu.hat, mu.hat + phi.hat*(Xt[-N]-mu.hat))
    sigma.hatsq = sum( (Xt-Xt.hat)^2 )/N
    BIC.obj = N*log(sigma.hatsq)+ 3*log(N) #6 always there
  }
  else{
    tau.vec = loc.ind*(1:N) #convert binary to CPT location
    tau = tau.vec[tau.vec>0] #keep CPT locations only
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
  return(-BIC.obj)
}

xt=Z.sim



tim5 = Sys.time()
BICGA = GA::ga(type="binary", fitness = BinSearch.BIC,
               nBits = Ts, maxiter = 10000, run = 1000,
               popSize = 200, monitor = F)
tim6 = Sys.time()

BICsol = (BICGA@solution)[1, ] ##can be 2 sols since X[1] is free, take 1st sol
BICsol[1] = 0
tBIC = as.vector( which(BICsol == 1) )

tmp1$bestfit
tmp1$bestchrom[1:(tmp1$bestchrom[1]+2)]
tmp2$overbestfit
tmp2$overbestchrom
BICGA@fitnessValue
tBIC

tim2 - tim1
tim4 - tim3
tim6 - tim5
