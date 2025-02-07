# changepointGA
Modified genetic algorithm (GA) special for changepoint detection in time series.

## Installation
You can install the version of changepointGA from CRAN:

```{r}
install.packages("changepointGA")
```

or the development version from Github:

```{r}
# install.packages("devtools")
devtools::install_github("mli171/changepointGA")
```

## Changepoint Detection Only

### An example of using the GA 
```{r}
##### Stationary time series with autocorrelation
N = 1000
betaT = c(0.5, -0.5, 0.3) # intercept, B, D
period = 30
XMatT = cbind(rep(1, N), cos(2*pi*(1:N)/period), sin(2*pi*(1:N)/period))
colnames(XMatT) = c("intercept", "Bvalue", "DValue")
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(N*Cp.prop)

Xt = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)

TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)

GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxgen       = 100000,
  maxconv      = 1000,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)

tim1 = Sys.time()
tmp1 = GA(ObjFunc=BinSearch.BIC, N=N, GA_param=GA_param, Xt=Xt)
tim2 = Sys.time()
```

### An example of using the island-based GA 
```{r}
##### Stationary time series with autocorrelation
N = 1000
betaT = c(0.5) # intercept
XMatT = matrix(1, nrow=N, ncol=1)
colnames(XMatT) = "intercept"
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(N*Cp.prop)

Xt = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)


## No parallel computing
IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 100,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = FALSE, ###
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)

tim3 = Sys.time()
tmp2 = IslandGA(ObjFunc=BinSearch.BIC, N=N, IslandGA_param, Xt=Xt)
tim4 = Sys.time()


## Parallel computing
IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 100,
  option       = "cp",
  monitoring   = FALSE,
  parallel     = TRUE, ###
  nCore        = 10,
  tol          = 1e-5,
  seed         = NULL
)

tim5 = Sys.time()
tmp3 = IslandGA(BinSearch.BIC, N=N, IslandGA_param, Xt=Xt)
tim6 = Sys.time()

tim4 - tim3
tim6 - tim5

tmp2$overbestfit
tmp3$overbestfit

tmp2$overbestchrom
tmp3$overbestchrom
```

## Changepoint Detection + Model order selection

### An example of using the GA 
```{r}
N = 1000
betaT = c(0.5, -0.5, 0.3) # intercept, B, D
period = 30
XMatT = cbind(rep(1, N), cos(2*pi*(1:N)/period), sin(2*pi*(1:N)/period))
colnames(XMatT) = c("intercept", "Bvalue", "DValue")
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(N*Cp.prop)

Xt = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)

p.range = list(ar=c(0,2), ma=c(0,2))

GA_param = list(
  popsize      = 200,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxgen       = 10000,
  maxconv      = 1000,
  option       = "both",
  monitoring   = FALSE,
  parallel     = TRUE,
  nCore        = 10,
  tol          = 1e-5,
  seed         = NULL
)

tim1 = Sys.time()
tmp1 = GA(ObjFunc=ARIMA.BIC.Order, N=N, GA_param, p.range=p.range, XMat=XMatT, Xt=Xt)
tim2 = Sys.time()

tim2 - tim1
tmp1$overbestfit
tmp1$overbestchrom
```

### An example of using the island-based GA 
```{r}
N = 1000
betaT = c(0.5, -0.5, 0.3) # intercept, B, D
period = 30
XMatT = cbind(rep(1, N), cos(2*pi*(1:N)/period), sin(2*pi*(1:N)/period))
colnames(XMatT) = c("intercept", "Bvalue", "DValue")
sigmaT = 1
phiT = c(0.5, -0.5)
thetaT = c(0.8)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(N*Cp.prop)

myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)

p.range = list(ar=c(0,2), ma=c(0,2))

IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 50,
  option       = "both",
  monitoring   = FALSE,
  parallel     = FALSE,
  nCore        = NULL,
  tol          = 1e-5,
  seed         = NULL
)

tim3 = Sys.time()
tmp2 = IslandGA(ObjFunc=ARIMA.BIC.Order, N=N, IslandGA_param, p.range=p.range, XMat=XMatT, Xt=Xt)
tim4 = Sys.time()


IslandGA_param = list(
  popsize      = 40,
  Islandsize   = 5,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 10/N,
  minDist      = 1,
  mmax         = N/2 - 1,
  lmax         = 2 + N/2 - 1,
  maxMig       = 500,
  maxgen       = 100,
  maxconv      = 50,
  option       = "both",
  monitoring   = FALSE,
  parallel     = TRUE,
  nCore        = 5,
  tol          = 1e-5,
  seed         = NULL
)

tim5 = Sys.time()
tmp3 = IslandGA(ObjFunc=ARIMA.BIC.Order, N=N, IslandGA_param, p.range=p.range, XMat=XMatT, Xt=Xt)
tim6 = Sys.time()


tim4 - tim3
tim6 - tim5

tmp2$overbestfit
tmp3$overbestfit

tmp2$overbestchrom
tmp3$overbestchrom
```
