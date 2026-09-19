# changepointGA

Authors: Mo Li (mo.li@louisiana.edu), QiQi Lu (qlu2@vcu.edu)

---

## Overview

Detecting changepoints in a time series of length $N$ entails evaluating up to 
$2^{N-1}$ possible changepoint models, making exhaustive enumeration 
computationally infeasible. Genetic algorithms (GAs) provide a stochastic 
way to identify the structural changes: a population of candidate models 
evolves via selection, crossover, and mutation operators until it converges 
on one changepoint model that balances the goodness-of-fit with parsimony. 

The R package `changepointGA` represents each candidate model using an
integer-valued chromosome and provides three GA-based search strategies:

- `cptga()`: a basic single-population genetic algorithm;
- `cptgaisl()`: an island model genetic algorithm (IMGA);
- `cptgascisl()`: a structured and consensus-guided island model genetic
  algorithm (SC-IMGA).
  
The SC-IMGA extends the original island model algorithm through a structured
birth-death-relocate mutation operator, cross-island consensus guidance,
consensus-guided mutation candidate selection, and objective-based local
refinement. These components are designed to improve the search over
variable-dimensional changepoint configurations while retaining the general
objective-function framework of `changepointGA`.

Parallel computing is supported on multi-core hardware. Users may supply
custom objective functions or genetic operators, and the package supports
both changepoint detection alone and simultaneous changepoint detection and
integer-valued model-order selection.

The original `cptga()` and `cptgaisl()` algorithms have been evaluated through
extensive simulations and applications to array-based comparative genomic
hybridization data and a century-long temperature series. The SC-IMGA provides
an additional search strategy for more challenging multiple-changepoint
optimization problems.


## Installation
You can install the version of changepointGA from CRAN:

```r
install.packages("changepointGA")
```

or the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("mli171/changepointGA", build_vignettes = FALSE, force = TRUE)
```

## Changepoint Detection Only

Call the library.

```r
library(changepointGA)
```

### An example of using `cptga()`
```r
##### Stationary time series with autocorrelation
Ts = 1000
betaT = c(0.5, -0.5, 0.3) # intercept, B, D
period = 30
XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period))
colnames(XMatT) = c("intercept", "Bvalue", "DValue")
sigmaT = 1
phiT = c(0.5)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)

Xt = ts_sim(Ts=Ts, beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=NULL, 
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)

tim1 = Sys.time()
tmp1 = cptga(ObjFunc=arima_bic, N=Ts, XMat=XMatT, Xt=Xt)
tim2 = Sys.time()
summary(tmp1)
plot(tmp1, data=Xt)


tim2 - tim1
```

### An example of using `cptgaisl()`
```r
##### Stationary time series with autocorrelation
Ts = 1000
betaT = c(0.5) # intercept
XMatT = matrix(1, nrow=Ts, ncol=1)
colnames(XMatT) = "intercept"
sigmaT = 1
phiT = c(0.5)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)

Xt = ts_sim(Ts=Ts, beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=NULL, 
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)


## No parallel computing
tim3 = Sys.time()
tmp2 = cptgaisl(ObjFunc=arima_bic, N=Ts, XMat=XMatT, Xt=Xt)
tim4 = Sys.time()
summary(tmp2)
plot(tmp2, data=Xt)


## Parallel computing
tim5 = Sys.time()
tmp3 = cptgaisl(ObjFunc=arima_bic, N=Ts, parallel=TRUE, nCore=5, XMat=XMatT, Xt=Xt)
tim6 = Sys.time()
summary(tmp3)
plot(tmp3, data=Xt)

tim4 - tim3
tim6 - tim5
```

### An example of using `cptgascisl()`

```r
##### Stationary time series with autocorrelation
Ts = 1000
betaT = c(0.5) # intercept
XMatT = matrix(1, nrow=Ts, ncol=1)
colnames(XMatT) = "intercept"
sigmaT = 1
phiT = c(0.5)
DeltaT = c(2, -2)
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)

Xt = ts_sim(Ts=Ts, beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=NULL,
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)


## Structured and consensus-guided island model GA
tim7 = Sys.time()
tmp4 = cptgascisl(ObjFunc=arima_bic, N=Ts, XMat=XMatT, Xt=Xt)
tim8 = Sys.time()

summary(tmp4)
plot(tmp4, data=Xt)

tim8 - tim7
```

## Changepoint Detection + Model order selection

### An example of using `cptga()`
```r
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

Xt = ts_sim(Ts=Ts, beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)

prange = list(ar=c(0,2), ma=c(0,2))

tim1 = Sys.time()
tmp1 = cptga(ObjFunc=arima_bic_order_pq, N=Ts, prange=prange, option="both", 
             XMat=XMatT, Xt=Xt)
tim2 = Sys.time()
summary(tmp1)
plot(tmp1, data=Xt)

tim2 - tim1
```

### An example of using `cptgaisl()`
```r
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

Xt = ts_sim(Ts=Ts, beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, 
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)

prange = list(ar=c(0,2), ma=c(0,2))

tim3 = Sys.time()
tmp2 = cptgaisl(ObjFunc=arima_bic_order_pq, N=Ts, prange=prange, option="both", 
                XMat=XMatT, Xt=Xt)
tim4 = Sys.time()
summary(tmp2)
plot(tmp2, data=Xt)


tim5 = Sys.time()
tmp3 = cptgaisl(ObjFunc=arima_bic_order_pq, N=Ts, prange=prange, option="both", 
                parallel=TRUE, nCore=5, XMat=XMatT, Xt=Xt)
tim6 = Sys.time()
summary(tmp3)
plot(tmp3, data=Xt)

tim4 - tim3
tim6 - tim5
```

### An example of using `cptgascisl()`

```r
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

Xt = ts_sim(Ts=Ts, beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT,
            Delta=DeltaT, CpLoc=CpLocT, seed=1234)

prange = list(ar=c(0,2), ma=c(0,2))

tim7 = Sys.time()
tmp4 = cptgascisl(ObjFunc=arima_bic_order_pq, N=Ts, prange=prange,
                  option="both", XMat=XMatT, Xt=Xt)
tim8 = Sys.time()

summary(tmp4)
plot(tmp4, data=Xt)

tim8 - tim7
```

## Code style

Before pushing changes, please run 

```r
styler::style_pkg()
```

to ensure your code follows the tidyverse style guide.

## Citation

If you use `changepointGA` in your research, please cite the article corresponding to the method used:

- For `cptga()` and `cptgaisl()`, please cite:

  Li, M., & Lu, Q. (2026). *changepointGA: An R package for Fast Changepoint Detection via Genetic Algorithm*. The R Journal.

- For `cptgascisl()`, please cite:

  Li, M. (2026). *Structured and Consensus-Guided Island Model Genetic Algorithm for Multiple Changepoint Detection*. Manuscript prepared for publication.

