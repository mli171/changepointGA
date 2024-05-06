#' Time series simulation with changepoint effects
#'
#' This is a function to simulate time series \eqn{Z_{t}, t=1,\ldots,T_{s}} from a class of models,
#'  \deqn{Z_{t}=\mu_{t}+e_{t}.}
#' \itemize{
#'  \item{Stationary time series without autocorrelation} \cr \eqn{\mu_{t}} is a constant and \eqn{e_{t}}'s are independent and identically distributed as \eqn{N(0,\sigma)}.
#'  \item{Stationary time series with autocorrelation} \cr \eqn{\mu_{t}} is a constant
#'  and \eqn{e_{t}} follows an AR(1) process,
#'  \deqn{e_{t} = \phi_{t-1}e_{t-1}+\epsilon_{t},}
#' where \eqn{\epsilon_{t}}'s are independent and identically distributed as \eqn{N(0,\sigma)}.
#'  \item{Stationary with seasonality and autocorrelation} \cr \eqn{\mu_{t}} is not a constant,
#'  \deqn{\mu_{t} = ASin(\frac{2\pi t}{S})+BSin(\frac{2\pi t}{S}),}
#'  and \eqn{e_{t}} follows an AR(1) process,
#'  \deqn{e_{t} = \phi_{t-1}e_{t-1}+\epsilon_{t},}
#' where \eqn{\epsilon_{t}}'s are independent and identically distributed as \eqn{N(0,\sigma)}.
#'  \item{Stationary with seasonality, trend, and autocorrelation} \cr \eqn{\mu_{t}} is not a constant,
#'  \deqn{\mu_{t} = ASin(\frac{2\pi t}{S})+BSin(\frac{2\pi t}{S}) + \alpha t,}
#'  and \eqn{e_{t}} follows an AR(1) process,
#'  \deqn{e_{t} = \phi_{t-1}e_{t-1}+\epsilon_{t},}
#' where \eqn{\epsilon_{t}}'s are independent and identically distributed as \eqn{N(0,\sigma)}.
#' }
#' The changepoint effects could be introduced through the \eqn{\mu_{t}} as
#'  \deqn{\mu_{t} = \Delta_{1}I_{t>\tau_{1}} + \ldots + \Delta_{m}I_{t>\tau_{m}},}
#' where \eqn{1\leq\tau_{1}<\ldots<\tau_{m}\leq T_{s}} are the changepoint locations and
#' \eqn{\Delta_{1},\ldots,\Delta_{m}} are the changepoint parameter that need to be estimated.
#' @param theta A parameter vector contains other mean function parameters without changepoint parameters.
#' @param XMat The covairates for time series mean function without changepoint indicators.
#' @param sigma The standard deviation for time series residuals \eqn{\epsilon_{t}}.
#' @param phi A single value for the autocorrelation parameter for AR(1) errors.
#' @param Delta The parameter vector contains the changepoint parameters for time series mean function.
#' @param CpLoc A vector contains the changepoint locations range from \eqn{1\leq\tau\leq T_{s}}.
#' @param seed The random seed for simulation reproducibility.
#' @return Returns the simulated time series with attributes:
#' \item{\code{Z}}{The simulated time series.}
#' \item{Attributes}{
#' \itemize{
#'  \item{\code{DesignX}} The covariates include all changepoint indicators.
#'  \item{\code{mu}} A vector includes the mean values for simulated time series sequences.
#'  \item{\code{theta}} The true parameter vector (including changepoint effects).
#'  \item{\code{CpLoc}} The true changepoint locations where we introduce the mean changing effects.
#'  \item{\code{seed}} The random seed used for this simulation.
#'  }
#' }
#'
#' @import Rcpp
#' @import stats
#' @useDynLib changepointGA
#' @export
#' @examples
#' Ts = 1000
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#' DeltaT = c(2, -2)
#' sigmaT = 1
#'
#' ##### M1: Stationary time series without autocorrelation
#' thetaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' myts = ts.sim(theta=thetaT, XMat=XMatT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
#' TsPlotCheck(myts, tau=CpLocT)
#'
#' ##### M2: Stationary time series with autocorrelation
#' phiT = 0.5
#' thetaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' myts = ts.sim(theta=thetaT, XMat=XMatT, phi=phiT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
#' TsPlotCheck(myts, tau=CpLocT)
#'
#' ##### M3: Stationary with seasonality and autocorrelation
#' thetaT = c(0.5, -0.5, 0.3) # intercept, B, D
#' period = 30
#' XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period))
#' colnames(XMatT) = c("intercept", "Bvalue", "DValue")
#' myts = ts.sim(theta=thetaT, XMat=XMatT, phi=phiT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
#' TsPlotCheck(myts, tau=CpLocT)
#'
#' ##### M4: Stationary with seasonality, trend, and autocorrelation
#' # scaled trend if large number of sample size
#' thetaT = c(0.5, -0.5, 0.3, 0.01) # intercept, B, D, alpha
#' period = 30
#' XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period), 1:Ts)
#' colnames(XMatT) = c("intercept", "Bvalue", "DValue", "trend")
#' myts = ts.sim(theta=thetaT, XMat=XMatT, phi=phiT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
#' TsPlotCheck(myts, tau=CpLocT)
ts.sim = function(theta, XMat, sigma, phi=NULL, Delta=NULL, CpLoc=NULL, seed=NULL){

  if(!is.null(seed)){set.seed(seed)}

  Ts = NROW(XMat)

  if(length(theta) != NCOL(XMat)){stop("Number of covariates and parameter not equal!")}

  if(is.null(Delta)){
    if(is.null(CpLoc)){
      warnings("\n No changepoint effects!\n")
      DesignX = XMat
      mu = DesignX%*%theta
    }else{
      stop("Changepoint effects invalid!")
    }
  }else{
    if(is.null(CpLoc)){
      stop("Changepoint effects invalid!")
    }else{
      if(any(CpLoc>Ts) | any(CpLoc<0)){stop("Changepoint effects invalid!")}
      if(length(Delta) != length(CpLoc)){stop("Changepoint effects invalid!")}
      # changepoints
      tmptau = unique(c(CpLoc, Ts+1))
      CpMat = matrix(0, nrow=Ts, ncol=length(tmptau)-1)
      for(i in 1:NCOL(CpMat)){CpMat[tmptau[i]:(tmptau[i+1]-1),i] = 1}
      DesignX = cbind(XMat, CpMat)
      theta = c(theta, Delta)
      mu = DesignX%*%theta
    }
  }

  if(is.null(phi)){
    et = rnorm(n=Ts, mean=0, sd=sigma) # independent
  }else{
    et = arima.sim(n=Ts, list(ar=phi), sd=sigma) # sd argument is for WN
  }

  Z = et + mu

  attr(Z, 'DesignX') = DesignX
  attr(Z, 'mu') = mu
  attr(Z, 'theta') = theta
  attr(Z, 'CpEff') = Delta
  attr(Z, 'CpLoc') = CpLoc
  attr(Z, 'seed') = seed

  return(Z)
}


#' Plot the simulated time series
#'
#' This is a function to plot the simulated time series with segmentation visualization
#' by provided changepoint locations.
#' @param Z The simulated time series.
#' @param tau The provided changepoint locations.
#' @param mu The provided meam values for each time \eqn{t}.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
#' @examples
#' ##### M1: Stationary time series without autocorrelation
#' Ts = 1000
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#' DeltaT = c(2, -2)
#' sigmaT = 1
#' thetaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' myts = ts.sim(theta=thetaT, XMat=XMatT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
#'
#' TsPlotCheck(myts, tau=CpLocT)
TsPlotCheck = function(Z, tau=NULL, mu=NULL){

  Ts = length(Z)

  plot(x=1:Ts, y=Z, type="l", xlab="Time", ylab="Z")

  m = length(tau)

  if(!is.null(tau)){
    abline(v=tau, lty="dashed", col="blue", lwd=2)
  }else{
    cat("\n ---------- No changepoint specified ----------\n")
  }

  if(is.null(mu)){
    cat("\n Need to calculate sample mean of each segment.\n")
    # calculate the segment mean
    tauclc = c(1, tau, Ts+1)
    seg.len = diff(tauclc)
    ff = rep(0:m, times=seg.len)
    Xseg = split(Z, ff)
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    for(i in 1:(m+1)){
      segments(x0=tauclc[i], y0=mu.seg[i], x1=tauclc[i+1], y1=mu.seg[i], col="red", lty="dashed", lwd=2)
    }
  }else{
    lines(x=1:Ts, y=mu, col="red", lty="dashed", lwd=2)
  }

}
