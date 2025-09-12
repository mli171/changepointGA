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
  fmls <- names(formals(fn))
  fmls <- setdiff(fmls, "...")  # don't forward everything via "..."
  dots[intersect(names(dots), fmls)]
}