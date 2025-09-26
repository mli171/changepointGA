# changepointGA
Detecting changepoints in a time series of length $N$ entails evaluating up to 
$2^{N-1}$ possible changepoint models, making exhaustive enumeration 
computationally infeasible. Genetic algorithms (GAs) provide a stochastic 
way to identify the structural changes: a population of candidate models 
evolves via selection, crossover, and mutation operators until it converges 
on one changepoint model that balances the goodness-of-fit with parsimony. 
The R package changepointGA encodes each candidate model as an integer 
chromosome vector and supports both the basic single-population model GA 
and the island model GA. Parallel computing is implemented on multi-core 
hardware to further accelerate computation. Users may supply custom fitness 
functions or genetic operators, while a user-friendly wrapper streamlines 
routine analyses. Extensive simulations demonstrate that our package runs 
significantly faster than binary-encoded GA alternatives. Additionally, 
this package can simultaneously locate changepoints and estimate their 
effects, as well as other model parameters and any integer-valued 
hyperparameters. Applications to array-based comparative genomic 
hybridization data and a century-long temperature series further 
highlight the package’s value in biological and climate research.

## Installation
You can install the version of changepointGA from CRAN:

```{r}
install.packages("changepointGA")
```

or the development version from Github:

```{r}
# install.packages("remotes")
remotes::install_github("mli171/changepointGA", build_vignettes = FALSE, force = TRUE)
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

tim1 = Sys.time()
tmp1 = cptga(ObjFunc=ARIMA.BIC, N=N, XMat=XMatT, Xt=Xt)
tim2 = Sys.time()
summary(tmp1)
plot(tmp1, data=Xt)


tim2 - tim1
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


## No parallel computing
tim3 = Sys.time()
tmp2 = cptgaisl(ObjFunc=ARIMA.BIC, N=N, XMat=XMaT, Xt=Xt)
tim4 = Sys.time()
summary(tmp2)
plot(tmp2, data=Xt)


## Parallel computing
tim5 = Sys.time()
tmp3 = cptgaisl(ObjFunc=ARIMA.BIC, N=N, parallel=TRUE, nCore=5, XMat=XMaT, Xt=Xt)
tim6 = Sys.time()
summary(tmp3)
plot(tmp3, data=Xt)

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

p.range = list(ar=c(0,2), ma=c(0,2))

tim1 = Sys.time()
tmp1 = cptga(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, option="both", 
             XMat=XMatT, Xt=Xt)
tim2 = Sys.time()
summary(tmp1)
plot(tmp1, data=Xt)

tim2 - tim1
```

### An example of using the island model GA 
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

p.range = list(ar=c(0,2), ma=c(0,2))

tim3 = Sys.time()
tmp2 = cptgaisl(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, option="both", 
                XMat=XMatT, Xt=Xt)
tim4 = Sys.time()
summary(tmp2)
plot(tmp2, data=Xt)


tim5 = Sys.time()
tmp3 = cptgaisl(ObjFunc=ARIMA.BIC.Order, N=N, p.range=p.range, option="both", 
                parallel=TRUE, nCore=5, XMat=XMatT, Xt=Xt)
tim6 = Sys.time()
summary(tmp3)
plot(tmp3, data=Xt)

tim4 - tim3
tim6 - tim5
```

# Code style

Before pushing changes, please run 

```r
styler::style_pkg()
```

to ensure your code follows the tidyverse style guide.
