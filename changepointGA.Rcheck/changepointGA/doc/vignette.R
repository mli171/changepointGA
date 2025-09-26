## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(changepointGA)

## -----------------------------------------------------------------------------
Ts <- 200
betaT <- c(0.5) # intercept
XMatT <- matrix(rep(1, Ts), ncol = 1)
colnames(XMatT) <- c("intercept")
sigmaT <- 1
phiT <- c(0.5)
thetaT <- c(0.8)
DeltaT <- c(2, -2)
CpLocT <- c(50, 150)

myts <- ts.sim(beta = betaT, XMat = XMatT, sigma = sigmaT, phi = phiT, theta = thetaT, Delta = DeltaT, CpLoc = CpLocT, seed = 1234)
str(myts)

## ----fig.align = "center", fig.height=4, fig.width=6--------------------------
plot(x = 1:Ts, y = myts, type = "l", xlab = "Time", ylab = "Y")
abline(v = CpLocT, lty = "dashed", col = "blue", lwd = 2)

## Segmentation sample mean calculation and plotting
m <- length(CpLocT)
tauclc <- c(1, CpLocT, Ts + 1)
seg.len <- diff(tauclc)
ff <- rep(0:m, times = seg.len)
Xseg <- split(myts, ff)
mu.seg <- unlist(lapply(Xseg, mean), use.names = F)
for (i in 1:(m + 1)) {
  segments(
    x0 = tauclc[i], y0 = mu.seg[i],
    x1 = tauclc[i + 1], y1 = mu.seg[i],
    col = "red", lty = "dashed", lwd = 2
  )
}

## -----------------------------------------------------------------------------
ARIMA.BIC.Order(chromosome = c(2, 1, 1, 50, 150, Ts + 1), plen = 2, XMat = XMatT, Xt = myts)

## -----------------------------------------------------------------------------
N <- Ts
prange <- list(ar = c(0, 2), ma = c(0, 2))

## -----------------------------------------------------------------------------
suggestions <- list(NULL, c(50), c(50, 150), c(50, 100, 150))

## -----------------------------------------------------------------------------
XMatEst <- matrix(1, nrow = N, ncol = 1)

## -----------------------------------------------------------------------------
res.changepointGA <- suppressWarnings(cptga(
  ObjFunc = ARIMA.BIC.Order,
  N = N,
  prange = prange,
  suggestions = suggestions,
  option = "both",
  XMat = XMatEst,
  Xt = myts
))
print(res.changepointGA)
summary(res.changepointGA)

## -----------------------------------------------------------------------------
tim1 <- Sys.time()
res.Island.changepointGA <- suppressWarnings(cptgaisl(
  ObjFunc = ARIMA.BIC.Order,
  N = N,
  prange = prange,
  popSize = 160,
  numIslands = 2,
  maxMig = 1000,
  maxgen = 50,
  maxconv = 20,
  option = "both",
  XMat = XMatEst,
  Xt = myts
))
tim2 <- Sys.time()
tim2 - tim1
print(res.Island.changepointGA)
summary(res.Island.changepointGA)
plot(res.Island.changepointGA, data = myts)

## -----------------------------------------------------------------------------
tim3 <- Sys.time()
res.Island.changepointGA <- suppressWarnings(cptgaisl(
  ObjFunc = ARIMA.BIC.Order,
  N = N,
  prange = prange,
  popSize = 160,
  numIslands = 2,
  maxMig = 1000,
  maxgen = 50,
  maxconv = 20,
  option = "both",
  parallel = TRUE,
  nCore = 2,
  XMat = XMatEst,
  Xt = myts
))

tim4 <- Sys.time()
tim4 - tim3
print(res.Island.changepointGA)
summary(res.Island.changepointGA)
plot(res.Island.changepointGA, data = myts)

## -----------------------------------------------------------------------------
true.tau <- c(50, 150)
tau.Island <- res.Island.changepointGA@overbestchrom
est.tau <- c(tau.Island[4:(4 + tau.Island[1] - 1)])
cptDist(tau1 = true.tau, tau2 = est.tau, N = N)

## -----------------------------------------------------------------------------
sessionInfo()

