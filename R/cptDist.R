#' Comparing multiple changepoint configurations by pairwise distance
#'
#' This function is to calculate the distance metric (Shi et al., 2022) between
#' two changepoint configurations, \eqn{C_{1}=\{\tau_{1}, \ldots, \tau_{m}\}}
#' and \eqn{C_{2}=\{\eta_{1}, \ldots, \eta_{k}\}}, where \eqn{\tau}'s are the
#' changepoint locations.
#'
#' The pairwise distance was proposed by Shi et al. (2022),
#' \deqn{d(C_{1}, C_{2})=|m-k|+ \min(A(C_{1}, C_{2})),}
#' where \eqn{m} is the number of changepoints in configuration \eqn{C_{1}} and
#' \eqn{k} is the number of changepoints in configuration \eqn{C_{2}}.
#' The term \eqn{\min(A(C_{1}, C_{2}))} reflects the cost of matching
#' changepoint locations between \eqn{C_{1}} and \eqn{C_{2}} and can be calculated
#' using the linear assignment method. Details can be found in Shi et al. (2022).
#' Note: if one configuration doesn't contain any changepoints (valued \code{NULL}),
#' the distance is defined as \eqn{|m-k|}. The function can also be used
#' to examine changepoint detection performance in simulation studies. Given the
#' true changepoint configuration, \eqn{C_{true}}, used in generating the time
#' series, the calculated distance between the estimated multiple changepoint
#' configuration, \eqn{C_{est}=\{\eta_{1}, \ldots, \eta_{k}\}}, and \eqn{C_{true}}
#' can be used to evaluate the performance of the detection algorithm.
#' @param tau1 A vector contains the changepoint locations for \eqn{C_{1}}
#' for comparison. A value \code{NULL} is required if there is no changepoint detected.
#' @param tau2 A vector contains the changepoint locations for  \eqn{C_{2}}
#' for comparison. A value \code{NULL} is required if there is no changepoint detected.
#' @param N The simulated time series sample size. Two changepoint configurations
#' should have the same \eqn{n} values.
#' @return
#' \item{dist}{The calculated distance.}
#' @references{
#'   Shi, X., Gallagher, C., Lund, R., & Killick, R. (2022).
#'   A comparison of single and multiple changepoint techniques for time series data.
#'   \emph{Computational Statistics & Data Analysis}, 170, 107433.
#' }
#' @useDynLib changepointGA
#' @export
#' @examples
#' N = 100
#'
#' # both tau1 and tau2 has detected changepoints
#' tau2 = c(25, 50, 75)
#' tau1 = c(20, 35, 70, 80, 90)
#' cptDist(tau1=tau1, tau2=tau2, N=N)
#'
#' # either tau1 or tau2 has zero detected changepoints
#' cptDist(tau1=tau1, tau2=NULL, N=N)
#' cptDist(tau1=NULL, tau2=tau2, N=N)
#'
#' # both tau1 and tau2 has zero detected changepoints
#' cptDist(tau1=NULL, tau2=NULL, N=N)
cptDist = function(tau1, tau2, N){

  m = length(tau1)
  k = length(tau2)

  if(is.null(tau1) & is.null(tau2)){
    dist.res = NA
    warning("Both configurations are NULL.")
  }else{
    if(is.null(tau1) | is.null(tau2)){
      ACC = 0
    }else{
      if(any(tau1 < 0 | tau1 > N)){stop("First changepoint configuration invalid.")}
      if(any(tau2 < 0 | tau2 > N)){stop("Second changepoint configuration invalid.")}

      costs = matrix(0, nrow=m, ncol=k)
      for(i in 1:m){
        for(j in 1:k){
          costs[i,j] = abs(tau1[i] - tau2[j])/N
        }
      }
      if(m > k){costs=t(costs)}
      y = clue::solve_LSAP(costs, maximum = FALSE)
      ACC = sum(costs[cbind(seq_along(y), y)])
    }
    dist.res = abs(m-k) + ACC
  }

  return(dist.res)
}

