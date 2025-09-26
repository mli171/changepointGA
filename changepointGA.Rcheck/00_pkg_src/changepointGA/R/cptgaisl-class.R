#' S4 Class for Island Model Genetic Algorithm-Based Changepoint Detection
#'
#' This class stores results and settings for the island model genetic algorithm (`cptgaisl`)
#' used in multiple changepoint detection and optional model order selection.
#'
#' @name cptgaisl-class
#' @title S4 Class Definition for `cptgaisl`
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
#'
#' @return An object of class \code{cptgaisl}
#' @seealso \code{\link{cptgaisl}}, \code{\link{cptgaisl-class}}, \code{\link{random_population}}, \code{\link{selection_linearrank}}, \code{\link{uniformcrossover}}, \code{\link{mutation}}.

setClassUnion("numericOrNULL", members = c("numeric", "NULL"))
setClassUnion("listOrNULL", members = c("list", "NULL"))

#' @rdname cptgaisl-class
#' @export
setClass("cptgaisl",
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
    convg = "numeric"
  ),
  package = "changepointGA"
)

setMethod("print", "cptgaisl", function(x, ...) str(x))


#' Print Summary for a `cptgaisl` Object
#'
#' Displays key information about the settings and results from a changepoint detection
#' procedure using the Island model Genetic Algorithm (IMGA) stored in a `cptgaisl` object. This includes
#' the algorithm configuration, population settings, optimization mode, and final
#' solution such as the number and location of changepoints and model parameters (if applicable).
#'
#' @param x An object of class \code{cptgaisl}, typically produced by a GA-based changepoint detection routine.
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
#' @seealso \code{\link{cptgaisl}}, \code{\link[=summary.cptgaisl]{summary}}, \code{\link{plot.cptgaisl}}
#'
#' @method print summary.cptgaisl
#' @export
#' @aliases print.summary.cptgaisl
print.summary.cptgaisl <- function(x, digits = getOption("digits"), max_display = 5, ...) {
  cat("###############################################\n")
  cat("#  Changepoint Detection via Island Model GA  #\n")
  cat("###############################################\n")

  cat("   Settings: \n")
  cat(paste("   Population size         = ", x@popSize, "\n"))
  cat(paste("   Number of Island        = ", x@numIslands, "\n"))
  cat(paste("   Island size             = ", x@Islandsize, "\n"))
  cat(paste("   Number of generations   = ", x@count, "\n"))
  cat(paste("   Number of migrations    = ", x@countMig, "\n"))
  cat(paste("   Crossover probability   = ", format(x@pcrossover, digits = digits), "\n"))
  cat(paste("   Mutation probability    = ", format(x@pmutation, digits = digits), "\n"))
  cat(paste("   Changepoint probability = ", format(x@pchangepoint, digits = digits), "\n"))
  cat(paste("   Task mode               = ", x@option, "\n"))
  cat(paste("   Parallel Usage          = ", x@parallel, "\n"))
  if (x@parallel) {
    cat(paste("   Number of thread      = ", x@nCore, "\n"))
  }
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
  cat("\n##### Island GA results ##### \n")
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

  #
  invisible()
}

#' Print method for objects of class \code{cptgaisl}
#'
#' @param object An object of class \code{cptgaisl}.
#' @param ... Additional arguments (ignored).
#' @rdname cptgaisl-class
#' @aliases summary.cptgaisl
#' @export
setMethod("summary", "cptgaisl", function(object, ...) {
  print.summary.cptgaisl(object, ...)
})


#' Plot Time Series with Detected Changepoints from a `cptgaisl` Object
#'
#' This function visualizes a univariate time series along with the changepoints
#' identified by a island model genetic algorithm, as represented by a `cptgaisl` object.
#' Vertical dashed lines mark changepoint locations, and segment means are shown as horizontal
#' dashed lines. The optimal fitness value and changepoint locations are
#' displayed as margin text.
#'
#' @param x An object of class \code{cptgaisl}, typically returned by an island model genetic algorithm based
#' changepoint detection procedure.
#' @param data A numeric vector representing the observed univariate time series.
#' @param main Optional main title for the plot.
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
#' @seealso \code{\link[=summary,cptgaisl-method]{summary}}, \code{\link{print.summary.cptgaisl}}
#'
#' @exportS3Method
plot.cptgaisl <- function(x,
                          data,
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
  } else {
    tau <- NULL
  }

  fit <- sprintf("%.3f", x@overbestfit)
  tau_vals <- if (!is.null(tau)) if (use_custom_X) XTickLab[tau] else tau else NULL
  changepoint_str <- if (!is.null(tau_vals)) {
    paste0("Changepoints: ", paste(tau_vals, collapse = ", "))
  } else {
    "Changepoint Locations: None"
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

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
    ...
  )
  if (!is.null(main)) {
    title(main = main, line = 3.5) # push title down a bit
  }

  if (!is.null(XTickPos) && use_custom_X) {
    axis(1, at = match(XTickPos, XTickLab), labels = XTickPos)
  } else {
    axis(1, at = pretty(plot_x), labels = pretty(plot_x))
  }

  if (!is.null(tau)) {
    cp_x <- if (use_custom_X) XTickLab[tau] else tau
    abline(v = cp_x, col = "blue", lty = "dashed", lwd = lwd)

    tau_full <- c(1, tau, Ts + 1)
    seg_len <- diff(tau_full)
    ff <- rep(0:m, times = seg_len)
    mu.seg <- tapply(data, ff, mean)

    for (i in seq_along(mu.seg)) {
      segments(
        x0 = if (use_custom_X) XTickLab[tau_full[i]] else tau_full[i],
        y0 = mu.seg[i],
        x1 = if (use_custom_X) XTickLab[tau_full[i + 1]] else tau_full[i + 1],
        y1 = mu.seg[i],
        col = "red", lty = "dashed", lwd = lwd
      )
    }
  }

  mtext(paste("Fitness:", fit), side = 3, line = 1.5, adj = 0, cex = cex.lab)
  mtext(changepoint_str, side = 3, line = 0.5, adj = 0, cex = cex.lab)
}
