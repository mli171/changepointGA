#' Calculating BIC for Multiple changepoint detection with AR and MA order selection
#'
#' The objective function for changepoint search in Autoregressive
#' moving average with model order selection via Bayesian Information Criterion
#' (BIC).
#'
#' @param chromosome A vector consists of the number of changepoints, the order
#' of AR component (refers to the number of lagged terms used to model the
#' current value of a time series), the order of MA component (refers to the
#' number of lagged error terms used to model the current value of a time
#' series), the changepoint locations, and a value of time series sample size
#' plus 1 ($N+1$) indicating the end of the chromosome.
#' @param plen The number of model order parameters that need to be selected.
#' If model order selection needs to be performed simultaneously with the
#' changepoint detection task, \code{plen} should be nonzero.
#' @param XMat A matrix contains the covariates, but not includes changepoint
#' effects, for time series regression.
#' @param Xt The simulated ARMA time series from \code{ts.sim} function.
#' @return The BIC value of the objective function.
#' @import stats
#' @importFrom utils tail
#' @useDynLib changepointGA
#' @export
#' @examples
#' Ts <- 1000
#' XMatT <- matrix(1, nrow = Ts, ncol = 1)
#' Xt <- ts_sim(
#'   Ts = Ts, beta = 0.5, XMat = XMatT, sigma = 1, phi = 0.5, theta = 0.8,
#'   Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
#' )
#'
#' # one chromosome representation
#' chromosome <- c(2, 1, 1, 250, 750, 1001)
#' arima_bic_order(chromosome, plen = 2, XMat = XMatT, Xt = Xt)
arima_bic_order <- function(chromosome, plen = 2, XMat, Xt) {
  m <- chromosome[1]
  porder <- chromosome[2:(plen + 1)]
  tau <- chromosome[(plen + 2):length(chromosome)]
  N <- length(Xt) # length of the series

  if (m == 0) {
    ## Case 1, Zero Changepoint
    DesignX <- XMat
    fit <- try(arima(Xt,
      order = c(porder[1], 0, porder[2]), xreg = DesignX, include.mean = F,
      optim.control = list(maxit = 50)
    ))
    if (inherits(fit, "try-error")) {
      bic_obj <- NA
    } else {
      bic_obj <- BIC(fit)
    }
  } else {
    tau <- tau[tau > 1 & tau < N + 1] # keep CPT locations only
    tmptau <- unique(c(tau, N))
    CpMat <- matrix(0, nrow = N, ncol = length(tmptau) - 1)
    for (i in seq_len(NCOL(CpMat))) {
      CpMat[(tmptau[i] + 1):tmptau[i + 1], i] <- 1
    }
    DesignX <- cbind(XMat, CpMat)
    fit <- try(arima(Xt,
      order = c(porder[1], 0, porder[2]), xreg = DesignX, include.mean = F,
      optim.control = list(maxit = 50)
    ))
    if (inherits(fit, "try-error")) {
      bic_obj <- NA
    } else {
      bic_obj <- BIC(fit)
    }
  }

  return(bic_obj)
}
