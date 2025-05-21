#' Time series simulation with changepoint effects
#'
#' This is a function to simulate time series with changepoint effects. See
#' details below.
#'
#' @param beta A parameter vector contains other mean function parameters without changepoint parameters.
#' @param XMat The covairates for time series mean function without changepoint indicators.
#' @param sigma The standard deviation for time series residuals \eqn{\epsilon_{t}}.
#' @param phi A vector for the autoregressive (AR) parameters for AR part.
#' @param theta A vector for the moving average (MA) parameters for MA part.
#' @param Delta The parameter vector contains the changepoint parameters for time series mean function.
#' @param CpLoc A vector contains the changepoint locations range from \eqn{1\leq\tau\leq T_{s}}.
#' @param seed The random seed for simulation reproducibility.
#' @details
#' The simulated time series \eqn{Z_{t}, t=1,\ldots,T_{s}} is from a class of models,
#'  \deqn{Z_{t}=\mu_{t}+e_{t}.}
#' \itemize{
#'  \item{Time series observations are IID} \cr \eqn{\mu_{t}} is a constant
#'  and \eqn{e_{t}}'s are independent and identically distributed as \eqn{N(0,\sigma)}.
#'  \item{ARMA(p,q) model with constant mean} \cr \eqn{\mu_{t}} is a constant,
#'  and \eqn{e_{t}} follows an ARMA(p,q) process.
#'  \item{ARMA(p,q) model with seasonality} \cr \eqn{\mu_{t}} is not a constant,
#'  \deqn{\mu_{t} = ASin(\frac{2\pi t}{S})+BSin(\frac{2\pi t}{S}),}
#'  and \eqn{e_{t}} follows an ARMA(p,q) process.
#'  \item{ARMA(p,q) model with seasonality and trend} \cr \eqn{\mu_{t}} is not a constant,
#'  \deqn{\mu_{t} = ASin(\frac{2\pi t}{S})+BSin(\frac{2\pi t}{S}) + \alpha t,}
#'  and \eqn{e_{t}} follows an ARMA(p,q) process.
#' The changepoint effects could be introduced through the \eqn{\mu_{t}} as
#'  \deqn{\mu_{t} = \Delta_{1}I_{t>\tau_{1}} + \ldots + \Delta_{m}I_{t>\tau_{m}},}
#' where \eqn{1\leq\tau_{1}<\ldots<\tau_{m}\leq T_{s}} are the changepoint locations and
#' \eqn{\Delta_{1},\ldots,\Delta_{m}} are the changepoint parameter that need to be estimated.
#' }
#' @return The simulated time series with attributes:
#' \item{\code{Z}}{The simulated time series.}
#' \item{Attributes}{
#' \itemize{
#'  \item{\code{DesignX}} The covariates include all changepoint indicators.
#'  \item{\code{mu}} A vector includes the mean values for simulated time series sequences.
#'  \item{\code{theta}} The true parameter vector (including changepoint effects).
#'  \item{\code{CpEff}} The true changepoint magnitudes.
#'  \item{\code{CpLoc}} The true changepoint locations where we introduce the mean changing effects.
#'  \item{\code{arEff}} A vector givies the AR coefficients.
#'  \item{\code{maEff}} A vector givies the MA coefficients.
#'  \item{\code{seed}} The random seed used for this simulation.
#'  }
#' }
#'
#' @import Rcpp
#' @import stats
#' @useDynLib changepointGA
#' @export
#' @examples
#' ##### M1: Time series observations are IID
#' Ts = 1000
#' betaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' sigmaT = 1
#' DeltaT = c(2, -2)
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#'
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT,
#'               Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#' TsPlotCheck(Y=myts, tau=CpLocT)
#'
#' ##### M2: ARMA(2,1) model with constant mean
#' Ts = 1000
#' betaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' sigmaT = 1
#' phiT = c(0.5, -0.5)
#' thetaT = c(0.8)
#' DeltaT = c(2, -2)
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#'
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT,
#'               phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#' TsPlotCheck(Y=myts, tau=CpLocT)
#'
#' ##### M3: ARMA(2,1) model with seasonality
#' Ts = 1000
#' betaT = c(0.5, -0.5, 0.3) # intercept, B, D
#' period = 30
#' XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period))
#' colnames(XMatT) = c("intercept", "Bvalue", "DValue")
#' sigmaT = 1
#' phiT = c(0.5, -0.5)
#' thetaT = c(0.8)
#' DeltaT = c(2, -2)
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#'
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT,
#'               phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#' TsPlotCheck(Y=myts, tau=CpLocT)
#'
#'
#' ##### M4: ARMA(2,1) model with seasonality and trend
#' # scaled trend if large number of sample size
#' Ts = 1000
#' betaT = c(0.5, -0.5, 0.3, 0.01) # intercept, B, D, alpha
#' period = 30
#' XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period), 1:Ts)
#' colnames(XMatT) = c("intercept", "Bvalue", "DValue", "trend")
#' sigmaT = 1
#' phiT = c(0.5, -0.5)
#' thetaT = c(0.8)
#' DeltaT = c(2, -2)
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#'
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT,
#'               phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#' TsPlotCheck(Y=myts, tau=CpLocT)
ts.sim = function(beta, XMat, sigma, phi=NULL, theta=NULL, Delta=NULL, CpLoc=NULL, seed=NULL){

  if(!is.null(seed)){set.seed(seed)}

  Ts = NROW(XMat)

  if(is.null(Delta)){
    if(is.null(CpLoc)){
      warnings("\n No changepoint effects!\n")
      DesignX = XMat
      mu = DesignX%*%beta
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
      beta = c(beta, Delta)
      mu = DesignX%*%beta
    }
  }

  if(is.null(phi) & is.null(theta)){
    et = rnorm(n=Ts, mean=0, sd=sigma) # independent
  }else{
    et = arima.sim(n=Ts, list(ar=phi, ma=theta), sd=sigma) # sd argument is for WN
  }

  Z = et + mu

  attr(Z, 'DesignX') = DesignX
  attr(Z, 'mu') = mu
  attr(Z, 'theta') = beta
  attr(Z, 'CpEff') = Delta
  attr(Z, 'CpLoc') = CpLoc
  attr(Z, 'arEff') = phi
  attr(Z, 'maEff') = theta
  attr(Z, 'seed')  = seed

  return(Z)
}


#' Plot the simulated time series
#'
#' This is a function to plot the simulated time series with segmentation visualization
#' by provided changepoint locations.
#' @param X The time series time index, which could be specified as years, months, days, or others.
#' The default value is NULL and the vector from 1 to the time series length will be applied.
#' @param Xat The values from \code{X} that will be used as the X axis tick marks.
#' @param Y The time series data.
#' @param tau The provided changepoint locations.
#' @param mu The provided meam values for each time \eqn{t}.
#' @param XLAB A descriptive label for X axis.
#' @param YLAB A descriptive label for Y axis.
#' @import Rcpp
#' @import stats
#' @import graphics
#' @useDynLib changepointGA
#' @export
#' @return No return value, called for side effects
#' @examples
#' Ts = 1000
#' betaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' sigmaT = 1
#' DeltaT = c(2, -2)
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#'
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#' TsPlotCheck(X=1:Ts, Xat=seq(from=1, to=Ts, length=10), Y=myts, tau=CpLocT)
TsPlotCheck = function(X=NULL, Xat=NULL, Y, tau=NULL, mu=NULL, XLAB=NULL, YLAB=NULL){

  Ts = length(Y)

  if(is.null(XLAB)){XLAB = "Time"}
  if(is.null(YLAB)){YLAB = "Y"}

  plot(x=1:Ts, y=Y, type="l", xlab=XLAB, ylab=YLAB, xaxt = "n")
  if (!is.null(Xat)){
    p = match(Xat, X)
    axis(1, at=p, labels=X[p])
  }else{
    axis(1, at=1:Ts, labels=1:Ts)
  }
  m = length(tau)

  if(!is.null(tau)){
    abline(v=tau, lty="dashed", col="blue", lwd=2)
  }else{
    message("\n ---------- No changepoint specified ----------\n")
  }

  if(is.null(mu)){
    # calculate the segment mean
    tauclc = c(1, tau, Ts+1)
    seg.len = diff(tauclc)
    ff = rep(0:m, times=seg.len)
    Xseg = split(Y, ff)
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    for(i in 1:(m+1)){
      segments(x0=tauclc[i], y0=mu.seg[i], x1=tauclc[i+1], y1=mu.seg[i], col="red", lty="dashed", lwd=2)
    }
  }else{
    lines(x=1:Ts, y=mu, col="red", lty="dashed", lwd=2)
  }

}
