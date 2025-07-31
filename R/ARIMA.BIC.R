#' Example function: Calculating BIC for AR(1) model
#'
#' The objective function for changepoint search in Autoregressive
#' moving average with model order selection via Bayesian Information Criterion
#' (BIC).
#'
#' @param chromosome The chromosome consists of the number of changepoints, the
#' order of AR part, the order of MA part, the changepoint locations, and a
#' value of time series length plus 1 (N+1) indicating the end of the chromosome.
#' @param plen The number of model order parameters that need to be selected.
#' If model order selection needs to be performed simultaneously with the
#' changepoint detection task, \code{plen} should be nonzero.
#' @param XMat A matrix contains the covariates, but not includes changepoint
#' effects, for time series regression.
#' @param Xt The simulated ARMA time series from \code{ts.sim} function.
#' @return The BIC value of the objective function.
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib changepointGA
#' @export
#' @examples
#' Ts <- 1000
#' betaT <- c(0.5) # intercept
#' XMatT <- matrix(1, nrow = Ts, ncol = 1)
#' colnames(XMatT) <- "intercept"
#' sigmaT <- 1
#' phiT <- c(0.5)
#' thetaT <- NULL
#' DeltaT <- c(2, -2)
#' Cp.prop <- c(1 / 4, 3 / 4)
#' CpLocT <- floor(Ts * Cp.prop)
#'
#' myts <- ts.sim(
#'   beta = betaT, XMat = XMatT, sigma = sigmaT, phi = phiT, theta = thetaT,
#'   Delta = DeltaT, CpLoc = CpLocT, seed = 1234
#' )
#'
#' # candidate changepoint configuration
#' chromosome <- c(2, 250, 750, 1001)
#' ARIMA.BIC(chromosome, XMat = XMatT, Xt = myts)
ARIMA.BIC <- function(chromosome, plen = 0, XMat, Xt) {
  ARIMA_BIC_changepointGA_rcpp(chromosome, XMat, Xt)
}
