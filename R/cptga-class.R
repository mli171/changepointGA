#' S4 Class for Genetic Algorithm-Based Changepoint Detection
#'
#' An object of class \code{cptga} stores results and configuration settings
#' for changepoint detection using a Genetic Algorithm (GA), optionally with
#' simultaneous model order selection. This class records GA control parameters,
#' intermediate population structures, and the optimal solution found.
#'
#' @name cptga-class
#' @title S4 Class Definition for `cptga`
#'
#' @slot call The matched call that created the object.
#' @slot N The sample size of the time series.
#' @slot p.range A list object. Default is \code{NULL}. If specified, it contains the
#' ranges for each model order parameter (integers). Required when \code{option = "both"}
#' is used for joint changepoint and model selection.
#' @slot popSize An integer representing the number of individuals in each GA population.
#' @slot pcrossover The probability that the crossover operator is applied to two chromosomes.
#' @slot pmutation The probability that the mutation operator is applied to a chromosome.
#' @slot pchangepoint The prior probability that a changepoint has occurred at each location.
#' @slot minDist The minimum allowed distance between two adjacent changepoints.
#' @slot mmax The maximum possible number of changepoints. Typically set based on time series length and \code{option}.
#' @slot lmax The maximum length of the chromosome. Typically set based on time series length and \code{option}.
#' @slot maxgen The maximum number of generations the GA is allowed to run.
#' @slot maxconv If the optimal fitness value does not improve over this many generations, GA stops.
#' @slot option A character string: either \code{"cp"} for changepoint detection only,
#' or \code{"both"} for changepoint and model order selection.
#' @slot monitoring Logical. If \code{TRUE}, prints intermediate GA progress.
#' @slot parallel Logical. If \code{TRUE}, enables parallel computation for fitness evaluation.
#' @slot nCore Integer or \code{NULL}. Number of cores used for parallel computation when \code{parallel = TRUE}.
#' @slot tol Numeric. Tolerance for determining GA convergence. Default is \code{1e-5}.
#' @slot seed An integer or \code{NULL}. Random seed for reproducibility.
#' @slot suggestions A list or \code{NULL}. Each element provides suggested changepoint locations
#' to guide initial population design and potentially accelerate convergence.
#' @slot population A matrix where each row represents an individual chromosome in the current population.
#' @slot fitness A numeric vector containing the fitness values of individuals in the current generation.
#' @slot overbestchrom A vector representing the best chromosome found over all generations.
#' @slot overbestfit A numeric scalar. The best (smallest) fitness value achieved.
#' @slot bestfit A numeric vector recording the best fitness value in each generation.
#' @slot count A numeric value indicating the number of generations the GA actually ran.
#' @slot convg A numeric vector representing convergence information. A value of \code{0} indicates the algorithm successful completion. A value of \code{1} indicates the the total number of generations exceeds the pre-specified \code{maxgen} limit.
#'
#' @return An object of class \code{cptga}.
#' @seealso \code{\link{GA_param}}, \code{\link{random_population}}, \code{\link{cptga-class}}



setClassUnion("numericOrNULL", members = c("numeric", "NULL"))
setClassUnion("listOrNULL", members = c("list", "NULL"))


#' @rdname cptga-class
#' @export
setClass(Class = "cptga", 
         representation(
           call = "language",
           N = "numeric",
           p.range = "listOrNULL",
           popSize = "numeric",
           pcrossover = "numeric",
           pmutation = "numeric",
           pchangepoint = "numeric",
           minDist = "numeric",
           mmax = "numericOrNULL",
           lmax = "numericOrNULL",
           maxgen = "numeric",
           maxconv = "numeric",
           option = "character",
           monitoring = "logical",
           parallel = "logical",
           nCore = "numericOrNULL",
           tol = "numeric",
           seed = "numericOrNULL",
           suggestions = "listOrNULL",
           population = "matrix",
           fitness = "vector",
           overbestchrom = "vector",
           overbestfit = "numeric",
           bestfit = "vector",
           count = "numeric",
           convg = "numeric"
         ),
         package = "changepointGA"
)

setMethod("print", "cptga", function(x, ...) str(x))

#' @exportS3Method
print.summary.cptga = function(x, digits=getOption("digits"), max_display=5, ...)
{
  cat("###############################################\n")
  cat("#         Changepoint Detection via GA        #\n")
  cat("###############################################\n")
  
  cat("   Settings: \n")
  cat(paste("   Population size         = ", x@popSize, "\n"))
  cat(paste("   Number of generations   = ", x@count, "\n"))
  cat(paste("   Crossover probability   = ", format(x@pcrossover, digits = digits), "\n"))
  cat(paste("   Mutation probability    = ", format(x@pmutation, digits = digits), "\n"))
  cat(paste("   Changepoint probability = ", format(x@pchangepoint, digits = digits), "\n"))
  cat(paste("   Task mode               = ", x@option, "\n"))
  cat(paste("   Parallel Usage          = ", x@parallel, "\n"))
  if(x@parallel){
    cat(paste("   Number of thread      = ", x@nCore, "\n"))
  }
  if(!is.null(x@suggestions)){ 
    cat("   Suggestions: \n")
    for (i in seq_along(x@suggestions)) {
      cat("    ", sprintf("[%d]:", i))
      cat(x@suggestions[[i]], sep = " ")
      if (i > max_display) {
        cat("\n     ...")
        break
      }
      cat("\n")
    }
  }
  cat("\n##### GA results ##### \n")
  cat(paste("   Optimal Fitness value =", format(x@overbestfit, digits = digits), "\n"))
  cat(paste("   Optimal Solution: \n")) 
  m.sol = x@overbestchrom[1]
  cat(paste("        Number of Changepoints = ", m.sol, "\n")) 
  if(x@option == "cp"){
    if(m.sol > 0){
      tau.sol = x@overbestchrom[2:(1+m.sol)]
      cat(paste("        Changepoints Locations = ", paste(tau.sol, collapse = " ")), "\n")
    }else{
      cat("        Changepoints Locations = No changepoint reached optimum \n")
    }
  }else if (x@option == "both"){
    n.hyperparam = length(x@p.range)
    name.hyperparam = names(x@p.range)
    if(is.null(name.hyperparam)){
      name.hyperparam = paste0("Hyper.param.", 1:n.hyperparam)
    }
    if(m.sol > 0){
      hyperparam.sol = x@overbestchrom[2:(1+n.hyperparam)]
      hyperparam.sol = paste0(name.hyperparam, " = ", hyperparam.sol)
      cat("        Model hyperparameters:\n")
      for(i in seq_along(name.hyperparam)){
        cat(paste("            ", hyperparam.sol[i]), "\n")
      }
      tau.sol = x@overbestchrom[(2+n.hyperparam):(1+n.hyperparam+m.sol)]
      cat(paste("        Changepoints Locations = ", paste(tau.sol, collapse = " ")), "\n")
    }else{
      hyperparam.sol = x@overbestchrom[2:(1+n.hyperparam)]
      hyperparam.sol = paste0(name.hyperparam, " = ", hyperparam.sol)
      cat("        Model hyperparameters:\n")
      for(i in seq_along(name.hyperparam)){
        cat(paste("            ", hyperparam.sol[i]), "\n")
      }
      cat("        Changepoints Locations = No changepoint reached optimum \n")
    }
  }
  
  #
  invisible()
}

#' @param object An object of class \code{cptga}.
#' @param ... Additional arguments passed to \code{print.summary.cptga}.
#'
#' @rdname cptga-class
#' @aliases summary,cptga-method
#' @export
setMethod("summary", "cptga", function(object, ...) {
  print.summary.cptga(object, ...)
})


# plot.cptga = function(X=NULL, Xat=NULL, Y, tau=NULL, mu=NULL, XLAB=NULL, YLAB=NULL){
#   
#   Ts = length(Y)
#   
#   if(is.null(XLAB)){XLAB = "Time"}
#   if(is.null(YLAB)){YLAB = "Y"}
#   
#   plot(x=1:Ts, y=Y, type="l", xlab=XLAB, ylab=YLAB, xaxt = "n")
#   if (!is.null(Xat)){
#     p = match(Xat, X)
#     axis(1, at=p, labels=X[p])
#   }else{
#     axis(1, at=1:Ts, labels=1:Ts)
#   }
#   m = length(tau)
#   
#   if(!is.null(tau)){
#     abline(v=tau, lty="dashed", col="blue", lwd=2)
#   }else{
#     message("\n ---------- No changepoint specified ----------\n")
#   }
#   
#   if(is.null(mu)){
#     # calculate the segment mean
#     tauclc = c(1, tau, Ts+1)
#     seg.len = diff(tauclc)
#     ff = rep(0:m, times=seg.len)
#     Xseg = split(Y, ff)
#     mu.seg = unlist(lapply(Xseg,mean), use.names=F)
#     for(i in 1:(m+1)){
#       segments(x0=tauclc[i], y0=mu.seg[i], x1=tauclc[i+1], y1=mu.seg[i], col="red", lty="dashed", lwd=2)
#     }
#   }else{
#     lines(x=1:Ts, y=mu, col="red", lty="dashed", lwd=2)
#   }
#   
# }
# 
# setMethod("plot", "cptga", plot.cptga)


