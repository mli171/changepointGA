
cptDist = function(tau1, tau2, n){

  m = length(tau1)
  k = length(tau2)

  if(is.null(tau1) & is.null(tau2)){
    dist = NA
    warning("Both configurations are NULL.")
  }else{
    if(is.null(tau1)){tau1 = rep(0,k)}
    tmpm = length(tau1)

    if(is.null(tau2)){tau2 = rep(0,m)}
    tmpk = length(tau2)

    if(any(tau1 < 0 | tau1 > n)){stop("First changepoint configuration invalid.")}
    if(any(tau2 < 0 | tau2 > n)){stop("Second changepoint configuration invalid.")}

    costs = matrix(0, nrow=tmpm, ncol=tmpk)
    for(i in 1:tmpm){
      for(j in 1:tmpk){
        costs[i,j] = abs(tau1[i] - tau2[j])/n
      }
    }
    y = clue::solve_LSAP(costs, maximum = FALSE)
    dist = abs(m-k) + sum(costs[cbind(seq_along(y), y)])
  }

  return(dist)
}

n = 100

tau1 = c(25, 50, 75)
tau2 = c(20, 35, 70, 80, 90)

cptDist(tau1=tau1, tau2=tau2, n=n)
cptDist(tau1=tau1, tau2=tau1, n=n)
cptDist(tau1=tau1, tau2=NULL, n=n)
cptDist(tau1=NULL, tau2=tau2, n=n)
cptDist(tau1=NULL, tau2=NULL, n=n)
