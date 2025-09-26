#' Internal: filter a dots-list to only args a function accepts
#'
#' @description
#' Returns the subset of `dots` whose names match the formal arguments of `fn`
#' (excluding `...`). Handy for routing `...` to different plugin functions.
#'
#' @param dots A named list (typically `list(...)`).
#' @param fn A function object or a character name of a function.
#'
#' @return A named list containing only the args that `fn` declares.
#' @keywords internal
#' @noRd
.filter_args <- function(dots, fn) {
  if (is.character(fn)) fn <- get(fn, mode = "function")
  stopifnot(is.function(fn))

  # guard: if dots unnamed or empty
  if (length(dots) == 0L || is.null(names(dots))) {
    return(list())
  }

  fmls <- names(formals(fn))
  fmls <- setdiff(fmls, "...")
  dots[intersect(names(dots), fmls)]
}

#' changepointGA: Genetic Algorithms for Changepoint Models
#'
#' Tools for fitting changepoint models with genetic algorithms and Armadillo-backed linear algebra.
#'
#' @keywords internal
#' @useDynLib changepointGA, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
