#' Example2
#'
#' The example objective function for changepoint search in Autoregressive
#' moving average (ARMA(2,1)) via Bayesian Information Criterion (BIC).
#'
#'
#' @param chromosome The chromosome consists of the number of changepoints, the
#' order of AR part, the order of MA part, the changepoint locations, and a
#' value of time series length plus 1 indicating the end of the chromosome.
#' @param plen Since we want to perform the model order selection and the changepoint
#' detection task, \code{plen} is nonzero and is set as the number of model
#' order parameter, which is 2.
#' @param XMat A matrix contains the covariates (not includes changepoint effects)
#' for time series regression.
#' @param Xt The simulated ARMA(2,1) time series from \code{ts.sim} function.
#' @return Returned the value of the obejctive function (i.e. BIC).
#' @import stats
#' @importFrom utils tail
#' @useDynLib changepointGA
#' @export
#' @examples
#' Ts = 1000
#' betaT = c(0.5) # intercept
#' XMatT = matrix(1, nrow=Ts, ncol=1)
#' colnames(XMatT) = "intercept"
#' sigmaT = 1
#' phiT = c(0.5)
#' thetaT = NULL
#' DeltaT = c(2, -2)
#' Cp.prop = c(1/4, 3/4)
#' CpLocT = floor(Ts*Cp.prop)
#'
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT,
#'               Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#'
#' # candidate changepoint configuration
#' chromosome = c(2, 2, 1, 250, 750, 1001)
#' ARIMASearch.BIC(chromosome, plen=2, XMat=XMatT, Xt=myts)
ARIMASearch.BIC = function(chromosome, plen=2, XMat, Xt){

  # tau here only contain the changepoint locations
  m = chromosome[1]
  p = chromosome[2:(plen+1)]
  tau = chromosome[(plen+2):length(chromosome)]

  N = length(Xt) #length of the series

  if(tail(tau, 1) != N+1 | length(tau) != m+1){
    cat("\n Current chromosome representation: \n", "\n", chromosome)
    stop("\n Error in chromosome representations!\n")
  }

  if(m == 0){
    ##Case 1, Zero Changepoint
    DesignX = XMat
    fit = arima(Xt, order = c(p[1],0,p[2]), xreg=DesignX, include.mean=F)
  }else{
    tau = tau[tau>1 & tau<N+1] #keep CPT locations only
    tmptau = unique(c(tau, N+1))
    CpMat = matrix(0, nrow=N, ncol=length(tmptau)-1)
    for(i in 1:NCOL(CpMat)){CpMat[tmptau[i]:(tmptau[i+1]-1),i] = 1}
    DesignX = cbind(XMat, CpMat)
    fit = arima(Xt, order=c(p[1],0,p[2]), xreg=DesignX, include.mean=F)
  }

  BIC.obj = -2*fit$loglik + (m+sum(p)+2)*log(N)

  return(BIC.obj)
}
