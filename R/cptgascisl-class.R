#' S4 Class for Structured and Consensus-Guided Island Model Genetic Algorithm-Based Changepoint Detection
#'
#' This class stores results and settings for the structured and consensus-guided island model genetic algorithm (`cptgascisl`)
#' used in multiple changepoint detection and optional model order selection.
#'
#' @name cptgascisl-class
#' @title S4 Class Definition for `cptgascisl`
#'
#' @slot call The matched call that created the object.
#' @slot N The sample size of the time series.
#' @slot prange A list or NULL. Ranges for each model order parameter when \code{option = "both"}.
#' @slot popSize Integer. It represents the total number of individuals in each generation, which equal to the number of islands multiplied by the size of each island (i.e., \code{popSize = numIslands × Islandsize}).
#' @slot numIslands Integer. The number of islands (sub-populations).
#' @slot Islandsize Numerical value. The number of individuals in each island.
#' @slot pcrossover Probability of crossover operation.
#' @slot pmutation Probability of mutation operation.
#' @slot pchangepoint Prior probability of changepoint occurrence.
#' @slot minDist Minimum distance between adjacent changepoints.
#' @slot mmax Maximum number of changepoints allowed.
#' @slot lmax Maximum length of chromosome.
#' @slot maxMig Maximum number of migrations allowed.
#' @slot maxgen Maximum number of generations per island before migration.
#' @slot maxconv Number of migrations with no improvement before stopping.
#' @slot option Either "cp" or "both".
#' @slot monitoring Logical. If TRUE, prints intermediate output.
#' @slot parallel Logical. Whether parallel computation is used.
#' @slot nCore Integer or NULL. Number of cores for parallelization.
#' @slot tol Tolerance threshold for fitness improvement.
#' @slot seed Integer or NULL. Seed for reproducibility.
#' @slot suggestions A list or NULL. Suggested changepoint configurations.
#' @slot Island A 3D array storing all individual chromosomes in current generation across islands. Dimensions are \code{lmax} × \code{Islandsize} × \code{numIslands}, representing chromosome length, individuals per island, and number of islands, respectively.
#' @slot IslandFit A matrix of fitness values in current generation with dimensions \code{Islandsize} × \code{numIslands}, where each column corresponds to one island's population.
#' @slot overbestchrom A vector. The best chromosome ever found.
#' @slot overbestfit Numeric. The best fitness score obtained.
#' @slot bestfit Numeric vector recording best fitness per migration.
#' @slot countMig Integer vector tracking number of migrations.
#' @slot count Integer vector tracking total generations.
#' @slot convg Integer vector for convergence diagnostics.
#' @slot consensus Logical. Whether cross-island consensus guidance is used.
#' @slot native_consensus Logical. Whether consensus guidance is used directly within the mutation operator.
#' @slot consensus_radius Nonnegative integer. Neighborhood radius used to construct the cross-island consensus score.
#' @slot consensus_lambda Numeric. Strength of consensus guidance.
#' @slot consensus_candidates Positive integer. Number of mutation proposals generated for consensus-guided candidate selection.
#' @slot local_refine Character. Objective refinement strategy: \code{"none"}, \code{"final"}, or \code{"periodic"}.
#' @slot local_every Positive integer. Number of migrations between periodic refinement steps.
#' @slot local_radius Nonnegative integer. Neighborhood radius used during objective refinement.
#' @slot local_max_passes Positive integer. Maximum number of objective refinement passes.
#' @slot n_eval_local Number of additional objective function evaluations performed during objective refinement.
#'
#' @return An object of class \code{cptgascisl}
#' @seealso \code{\link{cptgascisl}}, \code{\link{cptgascisl-class}}, \code{\link{random_population}}, \code{\link{selection_linear_rank}}, \code{\link{uniform_crossover}}, \code{\link{mutation_birth_death_relocate}}.

#' @rdname cptgascisl-class
#' @export
setClass("cptgascisl",
         representation(
           call = "language",
           N = "numeric",
           prange = "listOrNULL",
           popSize = "numeric",
           numIslands = "numeric",
           Islandsize = "numeric",
           pcrossover = "numeric",
           pmutation = "numeric",
           pchangepoint = "numeric",
           minDist = "numeric",
           mmax = "numeric",
           lmax = "numeric",
           maxMig = "numeric",
           maxgen = "numeric",
           maxconv = "numeric",
           option = "character",
           monitoring = "logical",
           parallel = "logical",
           nCore = "numericOrNULL",
           tol = "numeric",
           seed = "numericOrNULL",
           suggestions = "listOrNULL",
           Island = "array",
           IslandFit = "matrix",
           overbestchrom = "vector",
           overbestfit = "numeric",
           bestfit = "vector",
           countMig = "numeric",
           count = "numeric",
           convg = "numeric",
           consensus = "logical",
           native_consensus = "logical",
           consensus_radius = "numeric",
           consensus_lambda = "numeric",
           consensus_candidates = "numeric",
           local_refine = "character",
           local_every = "numeric",
           local_radius = "numeric",
           local_max_passes = "numeric",
           n_eval_local = "numeric"
         ),
         package = "changepointGA"
)


setMethod("print", "cptgascisl", function(x, ...) str(x))


#' Print Summary for a `cptgascisl` Object
#'
#' Displays key information about the settings and results from a changepoint detection
#' procedure using the Structured and Consensus-Guided Island Model Genetic Algorithm
#' (SC-IMGA) stored in a `cptgascisl` object. This includes the algorithm configuration,
#' population settings, consensus-guidance settings, objective-refinement settings,
#' optimization mode, and final solution such as the number and location of changepoints
#' and model parameters (if applicable).
#'
#' @param x An object of class \code{cptgascisl}, typically produced by a GA-based changepoint detection routine.
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
#' @seealso \code{\link{cptgascisl}}, \code{\link[=summary.cptgascisl]{summary}}, \code{\link{plot.cptgascisl}}
#'
#' @method print summary.cptgascisl
#' @export
#' @aliases print.summary.cptgascisl
print.summary.cptgascisl <- function(x, digits = getOption("digits"), max_display = 5, ...) {
  cat("#########################################################\n")
  cat("#  Changepoint Detection via Structured Consensus GA   #\n")
  cat("#########################################################\n")
  
  cat("   Settings: \n")
  cat(paste("   Population size         = ", x@popSize, "\n"))
  cat(paste("   Number of Island        = ", x@numIslands, "\n"))
  cat(paste("   Island size             = ", x@Islandsize, "\n"))
  cat(paste("   Number of generations   = ", x@count, "\n"))
  cat(paste("   Number of migrations    = ", x@countMig, "\n"))
  cat(paste("   Crossover probability   = ", format(x@pcrossover, digits = digits), "\n"))
  cat(paste("   Mutation probability    = ", format(x@pmutation, digits = digits), "\n"))
  cat(paste("   Changepoint probability = ", format(x@pchangepoint, digits = digits), "\n"))
  cat(paste("   minDist                 = ", x@minDist, "\n"))
  cat(paste("   Task mode               = ", x@option, "\n"))
  cat(paste("   Parallel Usage          = ", x@parallel, "\n"))
  if (x@parallel) {
    cat(paste("   Number of thread      = ", x@nCore, "\n"))
  }
  
  seed_print <- if (is.null(x@seed) || length(x@seed) == 0) {
    "NULL"
  } else {
    as.character(x@seed)
  }
  
  cat(paste("   Seed                    = ", seed_print, "\n"))
  
  if (!is.null(x@suggestions)) {
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
  
  cat("\n   Consensus settings: \n")
  cat(paste("   Consensus               = ", x@consensus, "\n"))
  cat(paste("   Native consensus        = ", x@native_consensus, "\n"))
  cat(paste("   Consensus radius        = ", x@consensus_radius, "\n"))
  cat(paste("   Consensus lambda        = ", format(x@consensus_lambda, digits = digits), "\n"))
  cat(paste("   Consensus candidates    = ", x@consensus_candidates, "\n"))
  
  cat("\n   Refinement settings: \n")
  cat(paste("   Local refinement        = ", x@local_refine, "\n"))
  cat(paste("   Local every             = ", x@local_every, "\n"))
  cat(paste("   Local radius            = ", x@local_radius, "\n"))
  cat(paste("   Local max passes        = ", x@local_max_passes, "\n"))
  cat(paste("   Local evaluations       = ", x@n_eval_local, "\n"))
  
  cat("\n##### Structured and Consensus-Guided Island Model GA results ##### \n")
  cat(paste("   Optimal Fitness value =", format(x@overbestfit, digits = digits), "\n"))
  cat(paste("   Optimal Solution: \n"))
  m.sol <- x@overbestchrom[1]
  cat(paste("        Number of Changepoints = ", m.sol, "\n"))
  
  if (x@option == "cp") {
    if (m.sol > 0) {
      tau.sol <- x@overbestchrom[2:(1 + m.sol)]
      cat(paste("        Changepoints Locations = ", paste(tau.sol, collapse = " ")), "\n")
    } else {
      cat("        Changepoints Locations = No changepoint reached optimum \n")
    }
  } else if (x@option == "both") {
    n.hyperparam <- length(x@prange)
    name.hyperparam <- names(x@prange)
    
    if (is.null(name.hyperparam)) {
      name.hyperparam <- paste0("Hyper.param.", 1:n.hyperparam)
    }
    
    if (m.sol > 0) {
      hyperparam.sol <- x@overbestchrom[2:(1 + n.hyperparam)]
      hyperparam.sol <- paste0(name.hyperparam, " = ", hyperparam.sol)
      cat("        Model hyperparameters:\n")
      for (i in seq_along(name.hyperparam)) {
        cat(paste("            ", hyperparam.sol[i]), "\n")
      }
      tau.sol <- x@overbestchrom[(2 + n.hyperparam):(1 + n.hyperparam + m.sol)]
      cat(paste("        Changepoints Locations = ", paste(tau.sol, collapse = " ")), "\n")
    } else {
      hyperparam.sol <- x@overbestchrom[2:(1 + n.hyperparam)]
      hyperparam.sol <- paste0(name.hyperparam, " = ", hyperparam.sol)
      cat("        Model hyperparameters:\n")
      for (i in seq_along(name.hyperparam)) {
        cat(paste("            ", hyperparam.sol[i]), "\n")
      }
      cat("        Changepoints Locations = No changepoint reached optimum \n")
    }
  }
  
  invisible()
}


#' Print method for objects of class \code{cptgascisl}
#'
#' @param object An object of class \code{cptgascisl}.
#' @param ... Additional arguments (ignored).
#' @rdname cptgascisl-class
#' @aliases summary.cptgascisl
#' @export
setMethod("summary", "cptgascisl", function(object, ...) {
  print.summary.cptgascisl(object, ...)
})


#' Plot Time Series with Detected Changepoints from a `cptgascisl` Object
#'
#' This function visualizes a univariate time series along with the changepoints
#' identified by a structured and consensus-guided island model genetic algorithm,
#' as represented by a `cptgascisl` object.
#' Vertical dashed lines mark changepoint locations, and segment means are shown as horizontal
#' dashed lines. The optimal fitness value and changepoint locations are
#' displayed as margin text.
#'
#' @param x An object of class \code{cptgascisl}, typically returned by a structured and consensus-guided island model genetic algorithm based
#' changepoint detection procedure.
#' @param data A numeric vector representing the observed univariate time series.
#' @param main Optional main title for the plot.
#' @param show_segmean Binary, whether to include the segments' means.
#' @param XTickLab Optional vector (e.g., numeric or date) for custom x-axis labels.
#'        Must be the same length as \code{data}.
#' @param XTickPos Optional vector specifying which elements of \code{XTickLab} to show as ticks.
#' @param XAxisLab Optional label for the x-axis. Default is \code{"Time"}.
#' @param YAxisLab Optional label for the y-axis. Default is \code{"Data"}.
#' @param cex.lab Text size for axis labels and margin text. Default is \code{1.3}.
#' @param cex.axis Text size for axis tick labels. Default is \code{1.3}.
#' @param cex.main Text size for the main title. Default is \code{1.3}.
#' @param lwd Line width for vertical and horizontal dashed lines. Default is \code{2}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @details
#' If \code{XTickLab} is supplied and matches the length of \code{data}, it is used for
#' the x-axis; otherwise, the default sequence \code{1:length(data)} is used.
#'
#' If the genetic algorithm was run with \code{option = "both"}, the function skips hyperparameters
#' in the chromosome when extracting changepoint positions.
#'
#' The plot displays vertical dashed lines at changepoint locations and horizontal dashed
#' lines for the mean of each segment. Fitness and changepoint summaries are shown above the plot.
#'
#' @return This function is called for its side effects and returns \code{NULL} invisibly.
#'
#' @seealso \code{\link[=summary,cptgascisl-method]{summary}}, \code{\link{print.summary.cptgascisl}}
#'
#' @exportS3Method
plot.cptgascisl <- function(x,
                            data,
                            show_segmean = TRUE,
                            main = NULL,
                            XTickLab = NULL,
                            XTickPos = NULL,
                            XAxisLab = "Time",
                            YAxisLab = "Data",
                            cex.lab = 1.3,
                            cex.axis = 1.3,
                            cex.main = 1.3,
                            lwd = 2, ...) {
  Ts <- length(data)
  use_custom_X <- !is.null(XTickLab) && length(XTickLab) == Ts
  plot_x <- if (use_custom_X) XTickLab else 1:Ts
  
  chrom <- x@overbestchrom
  m <- chrom[1]
  
  if (m > 0) {
    tau <- if (x@option == "both") {
      n.hyparam <- length(x@prange)
      chrom[(2 + n.hyparam):(1 + n.hyparam + m)]
    } else {
      chrom[2:(1 + m)]
    }
    tau <- sort(unique(tau))
  } else {
    tau <- integer(0)
  }
  
  fit <- sprintf("%.3f", x@overbestfit)
  
  starts <- c(1, tau + 1)
  ends   <- c(tau, Ts)
  
  ## changepoint labels
  if (length(tau) > 0) {
    tau_vals <- if (use_custom_X) XTickLab[tau] else tau
    changepoint_str <- paste0("Changepoints: ", paste(tau_vals, collapse = ", "))
  } else {
    changepoint_str <- "Changepoint Locations: None"
  }
  
  op <- par(c("mar", "cex.lab", "cex.axis", "cex.main"))
  on.exit(par(op), add = TRUE)
  
  par(
    mar = c(5, 5, 6, 2),
    cex.lab = cex.lab,
    cex.axis = cex.axis,
    cex.main = cex.main
  )
  
  plot(plot_x, data,
       type = "l",
       xlab = XAxisLab,
       ylab = YAxisLab,
       xaxt = "n",
       ...)
  
  if (!is.null(main)) {
    title(main = main, line = 3.5)
  }
  
  if (!is.null(XTickPos) && use_custom_X) {
    axis(1, at = XTickPos, labels = XTickPos)
  } else {
    axis(1, at = pretty(plot_x), labels = pretty(plot_x))
  }
  
  if (length(tau) > 0) {
    cp_x <- if (use_custom_X) XTickLab[tau] else tau
    abline(v = cp_x, col = "blue", lty = "dashed", lwd = lwd)
  }
  
  if (show_segmean) {
    mu.seg <- sapply(seq_along(starts), function(i) {
      mean(data[starts[i]:ends[i]])
    })
    
    for (i in seq_along(mu.seg)) {
      x0 <- if (use_custom_X) XTickLab[starts[i]] else starts[i]
      x1 <- if (use_custom_X) XTickLab[ends[i]]   else ends[i]
      
      segments(x0 = x0, y0 = mu.seg[i],
               x1 = x1, y1 = mu.seg[i],
               col = "red", lty = "dashed", lwd = lwd)
    }
  }
  
  mtext(paste("Fitness:", fit), side = 3, line = 1.5, adj = 0, cex = cex.lab)
  mtext(changepoint_str, side = 3, line = 0.5, adj = 0, cex = cex.lab)
}