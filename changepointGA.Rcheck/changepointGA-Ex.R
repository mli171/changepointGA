pkgname <- "changepointGA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "changepointGA-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('changepointGA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ARIMA.BIC.Order")
### * ARIMA.BIC.Order

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ARIMA.BIC.Order
### Title: Calculating BIC for Multiple changepoint detection with model
###   order selection
### Aliases: ARIMA.BIC.Order

### ** Examples

N <- 1000
XMatT <- matrix(1, nrow = N, ncol = 1)
Xt <- ts.sim(
  beta = 0.5, XMat = XMatT, sigma = 1, phi = 0.5, theta = 0.8,
  Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
)

# one chromosome representation
chromosome <- c(2, 1, 1, 250, 750, 1001)
ARIMA.BIC.Order(chromosome, plen = 2, XMat = XMatT, Xt = Xt)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ARIMA.BIC.Order", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ARIMA.BIC")
### * ARIMA.BIC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ARIMA.BIC
### Title: Example function: Calculating BIC for AR(1) model
### Aliases: ARIMA.BIC

### ** Examples

Ts <- 1000
betaT <- c(0.5) # intercept
XMatT <- matrix(1, nrow = Ts, ncol = 1)
colnames(XMatT) <- "intercept"
sigmaT <- 1
phiT <- c(0.5)
thetaT <- NULL
DeltaT <- c(2, -2)
Cp.prop <- c(1 / 4, 3 / 4)
CpLocT <- floor(Ts * Cp.prop)

myts <- ts.sim(
  beta = betaT, XMat = XMatT, sigma = sigmaT, phi = phiT, theta = thetaT,
  Delta = DeltaT, CpLoc = CpLocT, seed = 1234
)

# candidate changepoint configuration
chromosome <- c(2, 250, 750, 1001)
ARIMA.BIC(chromosome, XMat = XMatT, Xt = myts)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ARIMA.BIC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cptDist")
### * cptDist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cptDist
### Title: Comparing multiple changepoint configurations by pairwise
###   distance
### Aliases: cptDist

### ** Examples

N <- 100

# both tau1 and tau2 has detected changepoints
tau2 <- c(25, 50, 75)
tau1 <- c(20, 35, 70, 80, 90)
cptDist(tau1 = tau1, tau2 = tau2, N = N)

# either tau1 or tau2 has zero detected changepoints
cptDist(tau1 = tau1, tau2 = NULL, N = N)
cptDist(tau1 = NULL, tau2 = tau2, N = N)

# both tau1 and tau2 has zero detected changepoints
cptDist(tau1 = NULL, tau2 = NULL, N = N)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cptDist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cptga")
### * cptga

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cptga
### Title: Genetic algorithm
### Aliases: cptga

### ** Examples

## No test: 

N <- 1000
XMatT <- matrix(1, nrow = N, ncol = 1)
Xt <- ts.sim(
  beta = 0.5, XMat = XMatT, sigma = 1, phi = 0.5, theta = NULL,
  Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
)

## Multiple changepoint detection without model order selection

# without suggestions
GA.res <- cptga(ObjFunc = ARIMA.BIC, N = N, XMat = XMatT, Xt = Xt)
summary(GA.res)
plot(GA.res, data = Xt)

# with suggestions
suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
GA.res <- cptga(ObjFunc = ARIMA.BIC, N = N, suggestions = suggestions, XMat = XMatT, Xt = Xt)
summary(GA.res)
plot(GA.res, data = Xt)


## Multiple changepoint detection with model order selection

prange <- list(ar = c(0, 3), ma = c(0, 3))

# without suggestions
GA.res <- cptga(
  ObjFunc = ARIMA.BIC.Order, N = N, prange = prange,
  option = "both", XMat = XMatT, Xt = Xt
)
summary(GA.res)
plot(GA.res, data = Xt)

# with suggestions
suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
GA.res <- cptga(
  ObjFunc = ARIMA.BIC.Order, N = N, prange = prange,
  suggestions = suggestions, option = "both", XMat = XMatT, Xt = Xt
)
summary(GA.res)
plot(GA.res, data = Xt)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cptga", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cptgaisl")
### * cptgaisl

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cptgaisl
### Title: Island model based genetic algorithm
### Aliases: cptgaisl

### ** Examples

## No test: 

N <- 1000
XMatT <- matrix(1, nrow = N, ncol = 1)
Xt <- ts.sim(
  beta = 0.5, XMat = XMatT, sigma = 1, phi = 0.5, theta = NULL,
  Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
)

## Multiple changepoint detection without model order selection

# without suggestions
GAISL.res <- cptgaisl(ObjFunc = ARIMA.BIC, N = N, XMat = XMatT, Xt = Xt)
summary(GAISL.res)
plot(GAISL.res, data = Xt)

# with suggestions
suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
GAISL.res <- cptgaisl(ObjFunc = ARIMA.BIC, N = N, suggestions = suggestions, XMat = XMatT, Xt = Xt)
summary(GAISL.res)
plot(GAISL.res, data = Xt)


## Multiple changepoint detection with model order selection

prange <- list(ar = c(0, 3), ma = c(0, 3))

# without suggestions
GAISL.res <- cptgaisl(
  ObjFunc = ARIMA.BIC.Order, N = N, prange = prange,
  option = "both", XMat = XMatT, Xt = Xt
)
summary(GAISL.res)
plot(GAISL.res, data = Xt)

# with suggestions
suggestions <- list(NULL, 250, c(250, 500), c(250, 625), c(250, 500, 750))
GAISL.res <- cptgaisl(
  ObjFunc = ARIMA.BIC.Order, N = N, prange = prange,
  suggestions = suggestions, option = "both", XMat = XMatT, Xt = Xt
)
summary(GAISL.res)
plot(GAISL.res, data = Xt)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cptgaisl", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ts.sim")
### * ts.sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ts.sim
### Title: Time series simulation with changepoint effects
### Aliases: ts.sim

### ** Examples

##### M1: Time series observations are IID
Ts <- 1000
betaT <- c(0.5) # intercept
XMatT <- matrix(1, nrow = Ts, ncol = 1)
colnames(XMatT) <- "intercept"
sigmaT <- 1
DeltaT <- c(2, -2)
Cp.prop <- c(1 / 4, 3 / 4)
CpLocT <- floor(Ts * Cp.prop)

myts <- ts.sim(
  beta = betaT, XMat = XMatT, sigma = sigmaT,
  Delta = DeltaT, CpLoc = CpLocT, seed = 1234
)

##### M2: ARMA(2,1) model with constant mean
Ts <- 1000
betaT <- c(0.5) # intercept
XMatT <- matrix(1, nrow = Ts, ncol = 1)
colnames(XMatT) <- "intercept"
sigmaT <- 1
phiT <- c(0.5, -0.5)
thetaT <- c(0.8)
DeltaT <- c(2, -2)
Cp.prop <- c(1 / 4, 3 / 4)
CpLocT <- floor(Ts * Cp.prop)

myts <- ts.sim(
  beta = betaT, XMat = XMatT, sigma = sigmaT,
  phi = phiT, theta = thetaT, Delta = DeltaT, CpLoc = CpLocT, seed = 1234
)

##### M3: ARMA(2,1) model with seasonality
Ts <- 1000
betaT <- c(0.5, -0.5, 0.3) # intercept, B, D
period <- 30
XMatT <- cbind(rep(1, Ts), cos(2 * pi * (1:Ts) / period), sin(2 * pi * (1:Ts) / period))
colnames(XMatT) <- c("intercept", "Bvalue", "DValue")
sigmaT <- 1
phiT <- c(0.5, -0.5)
thetaT <- c(0.8)
DeltaT <- c(2, -2)
Cp.prop <- c(1 / 4, 3 / 4)
CpLocT <- floor(Ts * Cp.prop)

myts <- ts.sim(
  beta = betaT, XMat = XMatT, sigma = sigmaT,
  phi = phiT, theta = thetaT, Delta = DeltaT, CpLoc = CpLocT, seed = 1234
)


##### M4: ARMA(2,1) model with seasonality and trend
# scaled trend if large number of sample size
Ts <- 1000
betaT <- c(0.5, -0.5, 0.3, 0.01) # intercept, B, D, alpha
period <- 30
XMatT <- cbind(rep(1, Ts), cos(2 * pi * (1:Ts) / period), sin(2 * pi * (1:Ts) / period), 1:Ts)
colnames(XMatT) <- c("intercept", "Bvalue", "DValue", "trend")
sigmaT <- 1
phiT <- c(0.5, -0.5)
thetaT <- c(0.8)
DeltaT <- c(2, -2)
Cp.prop <- c(1 / 4, 3 / 4)
CpLocT <- floor(Ts * Cp.prop)

myts <- ts.sim(
  beta = betaT, XMat = XMatT, sigma = sigmaT,
  phi = phiT, theta = thetaT, Delta = DeltaT, CpLoc = CpLocT, seed = 1234
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ts.sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
