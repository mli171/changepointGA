# changepointGA
Modified genetic algorithm (GA) special for changepoint detection in time series.

## Installation
```{r}
devtools::install_github(repo = "mli171/changepointGA")
```

## Changepoint Detection Only

### An example of using the GA 
```{r}
##### Stationary time series with autocorrelation
Ts = 1000
betaT = c(0.5) # intercept
XMatT = matrix(1, nrow=Ts, ncol=1)
colnames(XMatT) = "intercept"
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)

myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
              Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(myts, tau=CpLocT)

GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/Ts,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxgen       = 100000,
  maxconv      = 1000,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
ga_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank_cpp",
                    crossover  = "offspring_uniformcrossover_cpp",
                    mutation   = "mutation")

tim1 = Sys.time()
tmp1 = GA(ObjFunc=BinSearch.BIC, n=Ts, GA_param, ga_operators, Xt=myts)
tim2 = Sys.time()
```

### An example of using the island-based GA 
```{r}
##### Stationary time series with autocorrelation
Ts = 1000
betaT = c(0.5) # intercept
XMatT = matrix(1, nrow=Ts, ncol=1)
colnames(XMatT) = "intercept"
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)

myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
              Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(myts, tau=CpLocT)


## No parallel computing
IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/Ts,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 100,
  option       = "cp",
  monitoring   = TRUE,
  parallel     = FALSE, ###
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
IslandGA_operators = list(population = "random_population_cpp",
                          selection  = "selection_linearrank_cpp",
                          crossover  = "offspring_uniformcrossover_cpp",
                          mutation   = "mutation")

tim3 = Sys.time()
tmp2 = IslandGA(ObjFunc=BinSearch.BIC, n=Ts, IslandGA_param, IslandGA_operators, Xt=myts)
tim4 = Sys.time()


## Parallel computing
IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/Ts,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 100,
  option       = "cp",
  monitoring   = TRUE,
  parallel     = TRUE, ###
  nCore        = 10,
  tol          = 1e-5,
  seed         = NULL
)
IslandGA_operators = list(population = "random_population_cpp",
                          selection  = "selection_linearrank_cpp",
                          crossover  = "offspring_uniformcrossover_cpp",
                          mutation   = "mutation")

tim5 = Sys.time()
tmp3 = IslandGA(BinSearch.BIC, n=Ts, IslandGA_param, IslandGA_operators, Xt=myts)
tim6 = Sys.time()

tim4 - tim3
tim6 - tim5

tmp2$bestfit
tmp3$bestfit
tmp2$bestchrom
tmp3$bestchrom
```

## Changepoint Detection + Model order selection

### An example of using the GA 
```{r}
Ts = 1000
betaT = c(0.5, -0.5, 0.3) # intercept, B, D
period = 30
XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period))
colnames(XMatT) = c("intercept", "Bvalue", "DValue")
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)

myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(myts, tau=CpLocT)

p.range = list(ar=c(0,2), ma=c(0,2))


GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/Ts,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxgen       = 10000,
  maxconv      = 1000,
  option       = "both",
  monitoring   = TRUE,
  parallel     = TRUE,
  nCore        = 10,
  tol          = 1e-5,
  seed         = NULL
)

ga_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank_cpp",
                    crossover  = "offspring_uniformcrossover_cpp",
                    mutation   = "mutation")

tim1 = Sys.time()
tmp1 = GA(ObjFunc=ARIMASearch.BIC, n=Ts, GA_param, ga_operators,
          p.range=p.range, XMat=XMatT, Xt=myts)
tim2 = Sys.time()

tim2 - tim1
tmp1$overbestfit
tmp1$overbestchrom
```

### An example of using the island-based GA 
```{r}
Ts = 1000
betaT = c(0.5, -0.5, 0.3) # intercept, B, D
period = 30
XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period))
colnames(XMatT) = c("intercept", "Bvalue", "DValue")
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)

myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
              Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(myts, tau=CpLocT)

p.range = list(ar=c(0,2), ma=c(0,2))

IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/Ts,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 50,
  option       = "both",
  monitoring   = TRUE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)
IslandGA_operators = list(population = "random_population_cpp",
                          selection  = "selection_linearrank_cpp",
                          crossover  = "offspring_uniformcrossover_cpp",
                          mutation   = "mutation")

tim3 = Sys.time()
tmp2 = IslandGA(ObjFunc=ARIMASearch.BIC, n=Ts, IslandGA_param, IslandGA_operators,
                p.range=p.range, XMat=XMatT, Xt=myts)
tim4 = Sys.time()


IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/Ts,
  minDist      = 2,
  mmax         = Ts/2 - 1,
  lmax         = 2 + Ts/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 50,
  option       = "both",
  monitoring   = TRUE,
  parallel     = TRUE,
  nCore        = 5,
  tol          = 1e-5,
  seed         = NULL
)
IslandGA_operators = list(population = "random_population_cpp",
                          selection  = "selection_linearrank_cpp",
                          crossover  = "offspring_uniformcrossover_cpp",
                          mutation   = "mutation")
                          
tim5 = Sys.time()
tmp3 = IslandGA(ObjFunc=ARIMASearch.BIC, n=Ts, IslandGA_param, IslandGA_operators,
                p.range=p.range, XMat=XMatT, Xt=myts)
tim6 = Sys.time()


tim4 - tim3
tim6 - tim5

tmp2$bestfit
tmp3$bestfit

tmp2$bestchrom
tmp3$bestchrom
```
