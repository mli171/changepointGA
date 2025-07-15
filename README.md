# changepointGA
Genetic algorithms (GA) are stochastic search techniques designed to address 
combinatorial optimization problems by mimicking the principles of natural 
selection and evolution. GAs have proven effective in both single and multiple 
changepoint analyses within time series data, where each chromosome encodes 
the hyperparameters, number, and locations of changepoints, along with the 
associated model parameters. Starting with a population of potential changepoint 
configurations, GAsutilize genetic operators—selection, crossover, and mutation—to 
evolve toward solutions with enhanced fitness. This paper presents the R package 
changepointGA, which encodes changepoint chromosomes in an integer format. 
Furthermore, changepointGA facilitates the dynamic and simultaneous estimation 
of changepoint models hyperparameters, changepoint configurations, and model 
parameters, leading to more robust and accurate analyses. Several simulation 
studies and real-world applications are discussed to illustrate the package 
capabilities.

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

Xt = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)

TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)


tim1 = Sys.time()
tmp1 = cptga(ObjFunc=ARIMA.BIC, N=N, XMat=XMaT, Xt=Xt)
tim2 = Sys.time()
tim2 - tim1
summary(tmp1)
plot(tmp1)
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

Xt = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)


## No parallel computing
tim3 = Sys.time()
tmp2 = cptgaisl(ObjFunc=ARIMA.BIC, N=N, XMat=XMaT, Xt=Xt)
tim4 = Sys.time()
summary(tmp2)
plot(tmp2)


## Parallel computing
tim5 = Sys.time()
tmp3 = cptgaisl(ObjFunc=ARIMA.BIC, N=N, parallel=TRUE, nCore=5, XMat=XMaT, Xt=Xt)
tim6 = Sys.time()
summary(tmp3)
plot(tmp3)

tim4 - tim3
tim6 - tim5
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

Xt = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)

p.range = list(ar=c(0,2), ma=c(0,2))

tim1 = Sys.time()
tmp1 = cptga(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, option="both", 
             XMat=XMatT, Xt=Xt)
tim2 = Sys.time()
summary(tmp1)
plot(tmp1)

tim2 - tim1
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

myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
              Delta=DeltaT, CpLoc=CpLocT, seed=1234)
TsPlotCheck(X=1:N, Xat=seq(from=1, to=N, length=10), Y=Xt, tau=CpLocT)

p.range = list(ar=c(0,2), ma=c(0,2))

tim3 = Sys.time()
tmp2 = cptgaisl(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, option="both", 
                XMat=XMatT, Xt=Xt)
tim4 = Sys.time()
summary(tmp2)
plot(tmp2)


tim5 = Sys.time()
tmp3 = cptgaisl(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, option="both", 
                parallel=TRUE, nCore=5, XMat=XMatT, Xt=Xt)
tim6 = Sys.time()
summary(tmp3)
plot(tmp3)

tim4 - tim3
tim6 - tim5
```
