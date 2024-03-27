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


# check the time series plot with separating by changepoint
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


Ts = 1000
Cp.prop = c(1/4, 3/4)
CpLocT = floor(Ts*Cp.prop)
DeltaT = c(2, -2)

sigmaT = 1

##### M1: Stationary time series without autocorrelation
thetaT = c(0.5) # intercept

XMatT = matrix(1, nrow=Ts, ncol=1)
colnames(XMatT) = "intercept"
myts = ts.sim(theta=thetaT, XMat=XMatT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
TsPlotCheck(myts, tau=CpLocT)

##### M2: Stationary time series with autocorrelation
phiT = 0.5
thetaT = c(0.5) # intercept

XMatT = matrix(1, nrow=Ts, ncol=1)
colnames(XMatT) = "intercept"
myts = ts.sim(theta=thetaT, XMat=XMatT, phi=phiT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
TsPlotCheck(myts, tau=CpLocT)

##### M3: Stationary with seasonality and autocorrelation
thetaT = c(0.5, -0.5, 0.3) # intercept, B, D
period = 30

XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period))
colnames(XMatT) = c("intercept", "Bvalue", "DValue")
myts = ts.sim(theta=thetaT, XMat=XMatT, phi=phiT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
TsPlotCheck(myts, tau=CpLocT)

##### M4: Stationary with seasonality, trend, and autocorrelation
# scaled trend if large number of sample size
thetaT = c(0.5, -0.5, 0.3, 0.01) # intercept, B, D, alpha
period = 30

XMatT = cbind(rep(1, Ts), cos(2*pi*(1:Ts)/period), sin(2*pi*(1:Ts)/period), 1:Ts)
colnames(XMatT) = c("intercept", "Bvalue", "DValue", "trend")
myts = ts.sim(theta=thetaT, XMat=XMatT, phi=phiT, sigma=sigmaT, Delta=DeltaT, CpLoc=CpLocT)
TsPlotCheck(myts, tau=CpLocT)
