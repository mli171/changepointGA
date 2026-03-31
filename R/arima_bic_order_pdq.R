#' Calculate BIC for changepoint detection with ARIMA order selection
#'
#' Objective function for simultaneous selection of changepoint locations and
#' ARIMA model orders using the Bayesian Information Criterion (BIC).
#'
#' Each chromosome encodes the number of changepoints, the ARIMA orders
#' \eqn{(p,d,q)}, and candidate changepoint locations. For a given chromosome,
#' the function constructs the corresponding changepoint design matrix, combines
#' it with the user-supplied regression design matrix, fits an ARIMA model with
#' regressors, and returns the BIC value.
#'
#' @param chromosome A numeric vector encoding the candidate model. The
#'   chromosome is assumed to have the form
#'   \code{c(m, p, d, q, tau1, tau2, ..., tauK)}, where \code{m} is the number
#'   of changepoints, \code{p}, \code{d}, and \code{q} are the ARIMA orders,
#'   and the remaining entries are candidate changepoint locations.
#' @param plen An integer giving the number of model-order entries in the
#'   chromosome. The default is \code{3}, corresponding to \eqn{(p,d,q)}.
#' @param XMat A matrix of regression covariates excluding changepoint
#'   indicators.
#' @param Xt A numeric vector containing the observed time series.
#'
#' @details
#' The first element of \code{chromosome} specifies the number of changepoints.
#' The next \code{plen} elements specify the ARIMA orders, which are rounded to
#' integers. Remaining elements are treated as candidate changepoint locations.
#'
#' If changepoints are present, the function constructs a block-indicator matrix
#' whose columns represent successive changepoint segments. This matrix is then
#' appended to \code{XMat} to form the regression design.
#'
#' Before fitting, columns with zero variance or non-finite values are removed.
#' Remaining regressors are mean-centered. The ARIMA model is then fitted using
#' \code{\link[stats]{arima}} with \code{include.mean = FALSE}, and the
#' corresponding BIC is returned.
#'
#' If model fitting fails, the function returns a large penalty value
#' (\code{1e10}) so that the candidate is disfavored during optimization.
#'
#' @return A numeric scalar giving the BIC of the fitted candidate model. If the
#'   ARIMA fit fails, the function returns \code{1e10}.
#'
#' @import stats
#' @useDynLib changepointGA
#' @export
#'
#' @examples
#' Ts <- 200
#' XMatT <- matrix(1, nrow = Ts, ncol = 1)
#' colnames(XMatT) <- "intercept"
#'
#' Xt <- ts_sim(
#'   Ts = Ts, beta = 0.5, XMat = XMatT, sigma = 1,
#'   phi = 0.5, theta = 0.3, d = 1,
#'   Delta = c(2, -2), CpLoc = c(50, 150), seed = 1234
#' )
#'
#' ## chromosome = c(m, p, d, q, tau1, tau2, ...)
#' chromosome <- c(2, 1, 1, 1, 50, 150, Ts + 1)
#'
#' arima_bic_order_pdq(
#'   chromosome = chromosome,
#'   plen = 3,
#'   XMat = XMatT,
#'   Xt = Xt
#' )
arima_bic_order_pdq <- function(chromosome, plen = 3, XMat, Xt) {
  m <- as.integer(round(chromosome[1]))
  porder <- as.integer(round(chromosome[2:(plen + 1)]))
  tau <- chromosome[(plen + 2):length(chromosome)]
  N <- length(Xt)
  
  if (m == 0) {
    DesignX <- XMat
  } else {
    tau <- tau[seq_len(min(m, length(tau)))]
    tau <- as.integer(round(tau))
    tau <- tau[tau >= 1 & tau < N]
    
    if (length(tau) == 0) {
      DesignX <- XMat
    } else {
      tmptau <- sort(unique(c(tau, N)))
      CpMat <- matrix(0, nrow = N, ncol = length(tmptau) - 1)
      for (i in seq_len(ncol(CpMat))) {
        CpMat[(tmptau[i] + 1):tmptau[i + 1], i] <- 1
      }
      DesignX <- cbind(XMat, CpMat)
    }
  }
  
  if (!is.null(DesignX)) {
    keep <- apply(as.matrix(DesignX), 2, function(z) sd(z) > 0)
    DesignX <- as.matrix(DesignX)[, keep, drop = FALSE]
    if (ncol(DesignX) == 0) DesignX <- NULL
  }
  
  fit <- try(
    arima(
      Xt,
      order = c(porder[1], porder[2], porder[3]),
      xreg = DesignX,
      include.mean = FALSE,
      optim.control = list(maxit = 50)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(1e10)
  }
  
  BIC(fit)
}