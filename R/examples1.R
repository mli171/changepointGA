#' Example1
#'
#' The example objective function for changepoint search in First order Gaussian
#' autoregressions (AR(1)) via Bayesian Information Criterion (BIC).
#'
#' The AR order here is fixed to be 1 and we only want to conduct changepoint
#' detection in this example.
#'
#' @param chromosome The chromosome consists of the number of changepoints, their
#' locations, and a value of time series length plus 1 indicating the end of the
#' chromosome.
#' @param plen Since we only want to perform the changepoint detection task,
#' \code{plen} has to be set as 0.
#' @param Xt The simulated AR(1) time series from \code{ts.sim} function.
#' @return Returned the value of the obejctive function (i.e. BIC).
#' @import stats
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
#' myts = ts.sim(beta=betaT, XMat=XMatT, sigma=sigmaT, phi=phiT, theta=thetaT, Delta=DeltaT, CpLoc=CpLocT, seed=1234)
#'
#'# candidate changepoint configuration
#' chromosome = c(2, 250, 750, 1001)
#' BinSearch.BIC(chromosome, Xt=myts)
BinSearch.BIC = function(chromosome, plen=0, Xt){

  # tau here only contain the changepoint locations
  m = chromosome[1]
  tau = chromosome[2:length(chromosome)]

  N = length(Xt) #length of the series

  if(tail(tau, 1) != N+1 | length(tau) != m+1){
    cat("\n Current chromosome representation: \n", "\n", chromosome)
    stop("\n Error in chromosome representations!\n")
  }

  if(m == 0){
    ##Case 1, Zero Changepoint
    mu.hat = mean(Xt)
    phi.hat = sum((Xt-mu.hat)[-N]*(Xt-mu.hat)[-1])/sum((Xt-mu.hat)[-1]^2)
    Xt.hat = c(mu.hat, mu.hat + phi.hat*(Xt[-N]-mu.hat))
    sigma.hatsq = sum( (Xt-Xt.hat)^2 )/N
    BIC.obj = N*log(sigma.hatsq)+ 3*log(N) #6 always there
  }else{
    tau = tau[tau>1 & tau<N+1] #keep CPT locations only
    tau.ext = c(1,tau,(N+1)) #include CPT boundary 1 and N+1

    ## Split Xt to regimes/segments to
    ## compute phi.hat and sigma.hat.sq
    seg.len = diff(tau.ext) #length of each segments
    ff = rep(0:m, times=seg.len) ##create factors for segmentation
    Xseg = split(Xt, ff) ##Segmentation list
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    mu.hat = rep(mu.seg, seg.len)
    phi.hat = sum((Xt-mu.hat)[-N]*(Xt-mu.hat)[-1])/sum((Xt-mu.hat)[-1]^2)
    Xt.hat = c(mu.hat[1], mu.hat[-1] + phi.hat*(Xt[-N]-mu.hat[-N]))
    sigma.hatsq = sum( (Xt-Xt.hat)^2 )/N
    BIC.obj = N*log(sigma.hatsq) + (2*m + 3)*log(N)
  }
  return(BIC.obj)
}
