#' Time series simulation with changepoint effects
#'
#' This is a function to simulate time series with changepoint effects. See
#' details below.
#' 
#' @param Ts A integer indicating simulated time series length (sample size).
#' @param beta A parameter vector contains other mean function parameters without changepoint parameters.
#' @param XMat The covairates for time series mean function without changepoint indicators.
#' @param sigma The standard deviation for time series residuals \eqn{\epsilon_{t}}.
#' @param phi A vector for the autoregressive (AR) parameters for AR part.
#' @param theta A vector for the moving average (MA) parameters for MA part.
#' @param d A nonnegative integer for the order of differencing in the ARIMA model. Default is 0, which reduces to an ARMA model.
#' @param Delta The parameter vector contains the changepoint parameters for time series mean function.
#' @param CpLoc A vector contains the changepoint locations range from \eqn{1\leq\tau\leq T_{s}}.
#' @param seed The random seed for simulation reproducibility.
#' @details
#' The simulated time series \eqn{Z_t, t = 1, \ldots, T_s} is generated from
#' \deqn{Z_t = \mu_t + e_t.}
#'
#' The error process \eqn{e_t} can be:
#' \itemize{
#'   \item IID Gaussian noise when \code{phi = NULL}, \code{theta = NULL},
#'     and \code{d = 0}.
#'   \item An ARMA(p,q) process when \code{d = 0} and at least one of
#'     \code{phi} or \code{theta} is specified.
#'   \item An ARIMA(p,d,q) process when \code{d > 0}.
#' }
#'
#' The mean structure \eqn{\mu_t} is specified through \code{XMat} and
#' \code{beta}, and changepoint effects can be introduced as
#' \deqn{\mu_t = X_t^\top \beta + \Delta_1 I(t \ge \tau_1) + \cdots +
#' \Delta_m I(t \ge \tau_m),}
#' where \eqn{1 \le \tau_1 < \cdots < \tau_m \le T_s} are changepoint
#' locations and \eqn{\Delta_1, \ldots, \Delta_m} are changepoint magnitudes.
#'
#' @return The simulated time series with attributes:
#' \itemize{
#'   \item \code{DesignX}: The covariates including changepoint indicators.
#'   \item \code{mu}: The mean vector for the simulated time series.
#'   \item \code{theta}: The true parameter vector including changepoint effects.
#'   \item \code{CpEff}: The true changepoint magnitudes.
#'   \item \code{CpLoc}: The true changepoint locations.
#'   \item \code{arEff}: The AR coefficients.
#'   \item \code{maEff}: The MA coefficients.
#'   \item \code{diffEff}: The differencing order.
#'   \item \code{seed}: The random seed used for simulation.
#' }
#'
#' @import Rcpp
#' @import stats
#' @useDynLib changepointGA
#' @export
#'
#' @examples
#' ##### M1: Time series observations are IID
#' Ts <- 1000
#' betaT <- c(0.5)
#' XMatT <- matrix(1, nrow = Ts, ncol = 1)
#' colnames(XMatT) <- "intercept"
#' sigmaT <- 1
#' DeltaT <- c(2, -2)
#' Cpprop <- c(1 / 4, 3 / 4)
#' CpLocT <- floor(Ts * Cpprop)
#'
#' myts <- ts_sim(
#'   Ts = Ts, beta = betaT, XMat = XMatT, sigma = sigmaT,
#'   Delta = DeltaT, CpLoc = CpLocT, seed = 1234
#' )
#'
#' ##### M2: ARMA(2,1) model with constant mean
#' Ts <- 1000
#' betaT <- c(0.5)
#' XMatT <- matrix(1, nrow = Ts, ncol = 1)
#' colnames(XMatT) <- "intercept"
#' sigmaT <- 1
#' phiT <- c(0.5, -0.5)
#' thetaT <- c(0.8)
#' DeltaT <- c(2, -2)
#' CpLocT <- floor(Ts * c(1 / 4, 3 / 4))
#'
#' myts <- ts_sim(
#'   Ts = Ts, beta = betaT, XMat = XMatT, sigma = sigmaT,
#'   phi = phiT, theta = thetaT, d = 0,
#'   Delta = DeltaT, CpLoc = CpLocT, seed = 1234
#' )
#'
#' ##### M3: ARIMA(1,1,1) model with changepoint effects
#' Ts <- 1000
#' betaT <- c(0.5)
#' XMatT <- matrix(1, nrow = Ts, ncol = 1)
#' colnames(XMatT) <- "intercept"
#' sigmaT <- 1
#' phiT <- c(0.5)
#' thetaT <- c(-0.3)
#' DeltaT <- c(2, -2)
#' CpLocT <- floor(Ts * c(1 / 4, 3 / 4))
#'
#' myts <- ts_sim(
#'   Ts = Ts, beta = betaT, XMat = XMatT, sigma = sigmaT,
#'   phi = phiT, theta = thetaT, d = 1,
#'   Delta = DeltaT, CpLoc = CpLocT, seed = 1234
#' )
ts_sim <- function(Ts, beta, XMat, sigma, phi = NULL, theta = NULL, d = 0,
                   Delta = NULL, CpLoc = NULL, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.null(Delta)) {
    if (is.null(CpLoc)) {
      warning("\n No changepoint effects!\n")
      DesignX <- XMat
      mu <- DesignX %*% beta
    } else {
      stop("Changepoint effects invalid!")
    }
  } else {
    if (is.null(CpLoc)) {
      stop("Changepoint effects invalid!")
    } else {
      if (any(CpLoc > Ts) | any(CpLoc < 0)) {
        stop("Changepoint effects invalid!")
      }
      if (length(Delta) != length(CpLoc)) {
        stop("Changepoint effects invalid!")
      }
      # changepoints
      tmptau <- unique(c(CpLoc, Ts + 1))
      CpMat <- matrix(0, nrow = Ts, ncol = length(tmptau) - 1)
      for (i in 1:NCOL(CpMat)) {
        CpMat[tmptau[i]:(tmptau[i + 1] - 1), i] <- 1
      }
      DesignX <- cbind(XMat, CpMat)
      beta <- c(beta, Delta)
      mu <- DesignX %*% beta
    }
  }
  
  if (is.null(phi) & is.null(theta) & d == 0) {
    et <- rnorm(n = Ts, mean = 0, sd = sigma) # independent
  } else {
    p <- if (is.null(phi)) 0 else length(phi)
    q <- if (is.null(theta)) 0 else length(theta)
    
    et <- arima.sim(
      n = Ts,
      model = list(order = c(p, d, q), ar = phi, ma = theta),
      sd = sigma
    )
    
    et <- as.numeric(et)
    if (length(et) > Ts) {
      et <- tail(et, Ts)
    }
  }
  
  Z <- et + mu
  
  attr(Z, "DesignX") <- DesignX
  attr(Z, "mu") <- mu
  attr(Z, "theta") <- beta
  attr(Z, "CpEff") <- Delta
  attr(Z, "CpLoc") <- CpLoc
  attr(Z, "arEff") <- phi
  attr(Z, "maEff") <- theta
  attr(Z, "diffEff") <- d
  attr(Z, "seed") <- seed
  
  return(Z)
}
