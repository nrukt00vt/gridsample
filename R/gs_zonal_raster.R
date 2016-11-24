globalVariables(names=c(".SD"))

#' @author Forrest R. Stevens, \email{forrest.stevens@louisville.edu}
#' @title Zonal statistics calculator
#' 
#' @description
#' This function calculates zonal statstics across a raster layer, for each polygon in a rasterized polygon layer.
#'
#' 
#' @param x Raster* layer. The layer that zonal statistics should be calculated from.
#' @param z Raster* layer. A rasterized version of the zonal layer.
#' @param stat character. Name of statistic used to calculate a value across each polygon in \code{z}. ex. "mean", "sum".
#' @param digits numeric. Number of significant digits in zonal statistic output.
#' @param na.rm logical. Defines whether to remove \code{NA}
#' @param ... Other variables
#'  
#' @return Vector of values representing the calculated statistic for each polygon, sorted by the order of polygons in the polygon layer.
#' 
#' 
#' 
#' @export

gs_zonal_raster <- function(x, z, stat = "mean", digits = 1, na.rm = TRUE, ...) {
  
	fun <- match.fun(stat)
	vals <- getValues(x)
	zones <- round(getValues(z), digits = digits)
	rDT <- data.table(vals, z=zones)
	setkey(rDT, z)
	rDT[, lapply(.SD,FUN=fun,na.rm=na.rm), by=z]
}
