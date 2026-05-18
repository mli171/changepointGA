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

myts <- ts_sim(Ts=Ts, beta = betaT, XMat = XMatT, sigma = sigmaT, phi = phiT, theta = thetaT, Delta = DeltaT, CpLoc = CpLocT, seed = 1234)
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
arima_bic_order_pq(chromosome = c(2, 1, 1, 50, 150, Ts + 1), plen = 2, XMat = XMatT, Xt = myts)

## -----------------------------------------------------------------------------
N <- Ts
prange <- list(ar = c(0, 2), ma = c(0, 2))

## -----------------------------------------------------------------------------
suggestions <- list(NULL, c(50), c(50, 150), c(50, 100, 150))

## -----------------------------------------------------------------------------
XMatEst <- matrix(1, nrow = N, ncol = 1)

## -----------------------------------------------------------------------------
reschangepointGA <- suppressWarnings(cptga(
  ObjFunc = arima_bic_order_pq,
  N = N,
  prange = prange,
  suggestions = suggestions,
  option = "both",
  XMat = XMatEst,
  Xt = myts
))
print(reschangepointGA)
summary(reschangepointGA)

## -----------------------------------------------------------------------------
tim1 <- Sys.time()
resIslandchangepointGA <- suppressWarnings(cptgaisl(
  ObjFunc = arima_bic_order_pq,
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
print(resIslandchangepointGA)
summary(resIslandchangepointGA)
plot(resIslandchangepointGA, data = myts)

## -----------------------------------------------------------------------------
tim3 <- Sys.time()
resIslandchangepointGA <- suppressWarnings(cptgaisl(
  ObjFunc = arima_bic_order_pq,
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
print(resIslandchangepointGA)
summary(resIslandchangepointGA)
plot(resIslandchangepointGA, data = myts)

## -----------------------------------------------------------------------------
truetau <- c(50, 150)
tauIsland <- resIslandchangepointGA@overbestchrom
esttau <- c(tauIsland[4:(4 + tauIsland[1] - 1)])
cpt_dist(tau1 = truetau, tau2 = esttau, N = N)

## -----------------------------------------------------------------------------
sessionInfo()

