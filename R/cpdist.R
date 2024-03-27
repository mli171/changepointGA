n = 100

tau1 = c(25, 50, 75)
tau2 = c(20, 35, 70, 80, 90)

cptDist = function(tau1, tau2, n){

  m = length(tau1)
  if(is.null(tau2)){tau2 = rep(0,m)}
  k = length(tau2)

  if(any(tau1 < 0 | tau1 > n)){stop("First changepoint configuration invalid.")}
  if(any(tau2 < 0 | tau2 > n)){stop("Second changepoint configuration invalid.")}

  costs = matrix(0, nrow=m, ncol=k)
  for(i in 1:m){
    for(j in 1:k){
      costs[i,j] = abs(tau1[i] - tau2[j])/n
    }
  }
  y = clue::solve_LSAP(costs, maximum = FALSE)
  dist = abs(m-k) + sum(costs[cbind(seq_along(y), y)])

  return(dist)
}

tau1 = c(25, 50, 75)
tau2 = c(20, 35, 70, 80, 90)

cptDist(tau1=tau1, tau2=tau2, n=n)
cptDist(tau1=tau1, tau2=tau1, n=n)
cptDist(tau1=tau1, tau2=NULL, n=n)

