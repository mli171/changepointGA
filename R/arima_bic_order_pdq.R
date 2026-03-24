#' Calculating BIC for multiple changepoint detection with ARIMA order selection
#'
#' The objective function for changepoint search in autoregressive integrated
#' moving average models with simultaneous model order selection via the
#' Bayesian Information Criterion (BIC). This function is designed for the
#' case where the orders of the AR, differencing, and MA components are all
#' selected together with changepoint locations.
#'
#' @param chromosome A vector consisting of the number of changepoints, the
#' order of the AR component (refers to the number of lagged terms used to
#' model the current value of a time series), the differencing order, the
#' order of the MA component (refers to the number of lagged error terms used
#' to model the current value of a time series), the changepoint locations,
#' and a value of time series sample size plus 1 (\eqn{N+1}) indicating the
#' end of the chromosome. The chromosome is represented as
#' \code{c(m, p, d, q, tau1, tau2, ..., N + 1)}.
#' @param plen The number of model order parameters to be selected. For this
#' function, \code{plen = 3}, corresponding to the AR, differencing, and MA
#' orders.
#' @param XMat A matrix containing the covariates, excluding changepoint
#' effects, for time series regression.
#' @param Xt The simulated ARIMA time series generated from the
#' \code{ts_sim} function.
#'
#' @return The BIC value of the objective function. If model fitting fails,
#' \code{NA} is returned.
#'
#' @details
#' This function fits an ARIMA(\eqn{p,d,q}) model, where \eqn{p}, \eqn{d},
#' and \eqn{q} are extracted from the chromosome and coerced to integer values.
#' If changepoints are present, the corresponding indicator variables are
#' appended to the design matrix and included as regression effects in the
#' fitted model.
#'
#' @import stats
#' @useDynLib changepointGA
#' @export
#'
#' @examples
#' Ts <- 1000
#' XMatT <- matrix(1, nrow = Ts, ncol = 1)
#' Xt <- ts_sim(
#'   Ts = Ts, beta = 0.5, XMat = XMatT, sigma = 1,
#'   phi = 0.5, theta = 0.8, d = 1,
#'   Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
#' )
#'
#' # one chromosome representation: c(m, p, d, q, tau1, tau2, N + 1)
#' chromosome <- c(2, 1, 1, 1, 250, 750, 1001)
#' arima_bic_order_pdq(chromosome, plen = 3, XMat = XMatT, Xt = Xt)
arima_bic_order_pdq <- function(chromosome, plen = 3, XMat, Xt) {
  m <- chromosome[1]
  porder <- chromosome[2:(plen + 1)]   # c(p, d, q)
  tau <- chromosome[(plen + 2):length(chromosome)]
  N <- length(Xt)
  
  porder <- as.integer(round(porder))
  
  if (m == 0) {
    DesignX <- XMat
    fit <- try(arima(Xt,
                     order = c(porder[1], porder[2], porder[3]),
                     xreg = DesignX, include.mean = FALSE,
                     optim.control = list(maxit = 50)
    ), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      bic_obj <- NA
    } else {
      bic_obj <- BIC(fit)
    }
  } else {
    tau <- tau[tau > 1 & tau < N + 1]
    tmptau <- unique(c(tau, N))
    CpMat <- matrix(0, nrow = N, ncol = length(tmptau) - 1)
    for (i in seq_len(NCOL(CpMat))) {
      CpMat[(tmptau[i] + 1):tmptau[i + 1], i] <- 1
    }
    DesignX <- cbind(XMat, CpMat)
    
    fit <- try(arima(Xt,
                     order = c(porder[1], porder[2], porder[3]),
                     xreg = DesignX, include.mean = FALSE,
                     optim.control = list(maxit = 50)
    ), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      bic_obj <- NA
    } else {
      bic_obj <- BIC(fit)
    }
  }
  
  return(bic_obj)
}