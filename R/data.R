#' Daily CO2 Efflux Measurements within a Temperate Forest
#'
#' A dataset containing daily CO2 efflux levels as a dependent variable and
#' temperature as an independent variables within a temperate forest
#' for a full year.
#'
#' @format ## `ex_data_1D`
#' A data frame with 365 rows and 3 columns:
#' \describe{
#'  \item{Time}{The day of the year}
#'  \item{Temp}{The temperature}
#'  \item{CO2}{The carbon dioxide efflux}
#' }
#'
#' @source <https://doi.org/10.1007/s11104-017-3506-4>
#' @examples
#' data(ex_data_1D)
"ex_data_1D"

#' Spatial Distribution of Soil CO2 Efflux for CONUS
#'
#' A dataset containing spatially distributed CO2 efflux levels as a dependent
#' variable and temperature as an independent variables within CONUS.
#'
#' @format ## `ex_data_2D`
#' A data frame with 903 rows and 4 columns:
#' \describe{
#'  \item{X}{The longitude}
#'  \item{Y}{The latitude}
#'  \item{rTemp}{The temperature}
#'  \item{rRS}{The carbon dioxide efflux}
#' }
#'
#' @source <https://doi.org/10.1111/gcb.15666>
#' @examples
#' data(ex_data_2D)
"ex_data_2D"
