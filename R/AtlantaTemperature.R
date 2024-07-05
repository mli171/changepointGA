#' Atlanta airport yearly average temperatures dataset
#'
#' The annual mean temperature data from 1879 to 2013 at Atlanta, Georgia’s
#' Hartsfield International Airport weather station.
#'
#' @docType data
#'
#' @usage data(AtlantaTemperature)
#'
#' @keywords datasets
#'
#' @references Lund, R. B., Beaulieu, C., Killick, R., Lu, Q., & Shi, X. (2023).
#' Good practices and common pitfalls in climate time series changepoint techniques:
#' A review. *Journal of Climate*, 36(23), 8041-8057.
#'
#' @examples
#' data(AtlantaTemperature)
#' \donttest{Xt = AtlantaTemperature[,2]
#'     N = length(Xt)
#'     XMat = matrix(1, nrow=N, ncol=1)}
"AtlantaTemperature"
