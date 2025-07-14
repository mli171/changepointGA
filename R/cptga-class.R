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


#' Print Summary for a `cptga` Object
#'
#' Displays key information about the settings and results from a changepoint detection
#' procedure using the Genetic Algorithm (GA) stored in a `cptga` object. This includes
#' the algorithm configuration, population settings, optimization mode, and final
#' solution such as the number and location of changepoints and model parameters (if applicable).
#'
#' @param x An object of class \code{cptga}, typically produced by a GA-based changepoint detection routine.
#' @param digits Number of digits to print for probabilities and fitness. Default taken from \code{getOption("digits")}.
#' @param max_display Maximum number of suggested solutions to display if suggestions are provided.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' When the GA is run in \code{option = "cp"} mode, only changepoint locations are shown.
#' If \code{option = "both"}, the output includes the selected model hyperparameters along
#' with changepoint locations.
#'
#' The function uses plain text output to print a formatted summary to the console. If
#' \code{x@suggestions} is provided, only up to \code{max_display} suggestions will be shown.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of printing to the console.
#'
#' @seealso \code{\link[=summary,cptga-method]{summary}}, \code{\link{cptga-class}}, \code{\link{plot.cptga}}
#'
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

#' Print method for objects of class \code{cptga}
#'
#' @param object An object of class \code{cptga}.
#' @param ... Additional arguments (ignored).
#'
#' @rdname cptga-class
#' @aliases summary, cptga-method
#' @export
setMethod("summary", "cptga", function(object, ...) {
  print.summary.cptga(object, ...)
})


#' Plot Time Series with Detected Changepoints from a `cptga` Object
#'
#' This function visualizes a univariate time series along with the changepoints
#' identified by a Genetic Algorithm, as represented by a `cptga` object. Vertical
#' dashed lines mark changepoint locations, and segment means are shown as horizontal
#' dashed lines. A legend summarizes the optimal fitness value and changepoint locations.
#'
#' @param x An object of class \code{cptga}, typically returned by a GA-based
#' changepoint detection procedure.
#' @param data A numeric vector representing the observed univariate time series.
#' @param ... Additional optional graphical arguments, including:
#' \itemize{
#'   \item \code{XTickLab}: Optional numeric or date vector for custom x-axis labels. Must be the same length as \code{data}.
#'   \item \code{XTickPos}: Optional vector of positions on the x-axis where tick marks should appear.
#'   \item \code{XAxisLab}, \code{YAxisLab}: Optional character strings for x-axis and y-axis labels.
#' }
#' @param main Optional main title for the plot.
#' @param LegendPos Character string indicating the position of the legend.
#'        Accepts standard legend positions such as \code{"topright"}, \code{"bottom"},
#'        \code{"topleft"}, etc. Default is \code{"topright"}.
#'
#' @details
#' If \code{XTickLab} is supplied and matches the length of \code{data}, it is used for the x-axis;
#' otherwise, the default \code{1:length(data)} is used.
#'
#' If the GA was run with \code{option = "both"}, the function adjusts the changepoint indices
#' to skip over model hyperparameters in the chromosome.
#'
#' The plot legend reports the optimal fitness value and changepoint locations (transformed
#' by \code{XTickLab} if applicable). The y-axis is extended slightly to make room for the legend.
#'
#' @return A time series plot is displayed. The function returns \code{NULL} invisibly.
#'
#' @seealso \code{\link[=summary,cptga-method]{summary}}, \code{\link{print.summary.cptga}}
#'
#' @exportS3Method
plot.cptga <- function(x, data, ..., main = NULL, LegendPos = "topright") {
  dots = list(...)
  
  XTickLab = dots$XTickLab
  XTickPos = dots$XTickPos
  XAxisLab = dots$XAxisLab %||% "Time"
  YAxisLab = dots$YAxisLab %||% "Y"
  
  Ts = length(data)
  use_custom_X = !is.null(XTickLab) && length(XTickLab) == Ts
  
  chrom = x@overbestchrom
  if (x@option == "cp") {
    m = chrom[1]
    tau = if (m > 0) chrom[2:(1 + m)] else NULL
  } else if (x@option == "both") {
    m = chrom[1]
    n.hyparam = length(x@p.range)
    tau = if (m > 0) chrom[(2 + n.hyparam):(1 + n.hyparam + m)] else NULL
  }
  fit = format(x@overbestfit, digits = 4)
  
  tau_vals = if (!is.null(tau)) if (use_custom_X) XTickLab[tau] else tau else NULL
  changepoint_str = if (!is.null(tau_vals)) {
    paste0("Changepoints: ", paste(tau_vals, collapse = ", "))
  } else {
    "Changepoints Locations: None"
  }
  
  legend_text = c(
    paste("Fitness:", fit),
    changepoint_str
  )
  
  plot_x = if (use_custom_X) XTickLab else 1:Ts
  
  y_range = range(data, na.rm = TRUE)
  y_buffer = 0.15 * diff(y_range)
  y_min = y_range[1]
  y_max = y_range[2] + y_buffer
  
  op = par(no.readonly = TRUE)
  on.exit(par(op))
  
  plot(plot_x, data,
       type = "l",
       xlab = XAxisLab,
       ylab = YAxisLab,
       xaxt = "n",
       main = main,
       ylim = c(y_min, y_max)
  )
  
  if (!is.null(XTickPos) && use_custom_X) {
    axis(1, at = match(XTickPos, XTickLab), labels = XTickPos)
  } else {
    axis(1, at = pretty(plot_x), labels = pretty(plot_x))
  }
  
  if (!is.null(tau)) {
    abline(v = if (use_custom_X) XTickLab[tau] else tau, col = "blue", lty = "dashed", lwd = 2)
  }
  
  if (!is.null(tau)) {
    tau_full = c(1, tau, Ts + 1)
    seg_len = diff(tau_full)
    ff = rep(0:m, times = seg_len)
    mu.seg = tapply(data, ff, mean)
    for (i in seq_along(mu.seg)) {
      x0 = tau_full[i]
      x1 = tau_full[i + 1]
      segments(
        x0 = if (use_custom_X) XTickLab[x0] else x0,
        y0 = mu.seg[i],
        x1 = if (use_custom_X) XTickLab[x1] else x1,
        y1 = mu.seg[i],
        col = "red", lty = "dashed", lwd = 2
      )
    }
  }
  
  legend(LegendPos,
         legend = legend_text,
         bty = "o", box.lwd = 1, box.col = "black")
}
