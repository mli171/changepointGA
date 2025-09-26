#' Calculating BIC for Multiple changepoint detection with model order selection
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
#' N <- 1000
#' XMatT <- matrix(1, nrow = N, ncol = 1)
#' Xt <- ts.sim(
#'   beta = 0.5, XMat = XMatT, sigma = 1, phi = 0.5, theta = 0.8,
#'   Delta = c(2, -2), CpLoc = c(250, 750), seed = 1234
#' )
#'
#' # one chromosome representation
#' chromosome <- c(2, 1, 1, 250, 750, 1001)
#' ARIMA.BIC.Order(chromosome, plen = 2, XMat = XMatT, Xt = Xt)
ARIMA.BIC.Order <- function(chromosome, plen = 2, XMat, Xt) {
  m <- chromosome[1]
  p.order <- chromosome[2:(plen + 1)]
  tau <- chromosome[(plen + 2):length(chromosome)]
  N <- length(Xt) # length of the series

  if (m == 0) {
    ## Case 1, Zero Changepoint
    DesignX <- XMat
    fit <- try(arima(Xt,
      order = c(p.order[1], 0, p.order[2]), xreg = DesignX, include.mean = F,
      optim.control = list(maxit = 50)
    ))
    if (inherits(fit, "try-error")) {
      BIC.obj <- NA
    } else {
      BIC.obj <- BIC(fit)
    }
  } else {
    tau <- tau[tau > 1 & tau < N + 1] # keep CPT locations only
    tmptau <- unique(c(tau, N))
    CpMat <- matrix(0, nrow = N, ncol = length(tmptau) - 1)
    for (i in 1:NCOL(CpMat)) {
      CpMat[(tmptau[i] + 1):tmptau[i + 1], i] <- 1
    }
    DesignX <- cbind(XMat, CpMat)
    fit <- try(arima(Xt,
      order = c(p.order[1], 0, p.order[2]), xreg = DesignX, include.mean = F,
      optim.control = list(maxit = 50)
    ))
    if (inherits(fit, "try-error")) {
      BIC.obj <- NA
    } else {
      BIC.obj <- BIC(fit)
    }
  }

  return(BIC.obj)
}
