##  Copyright (c) 2016, Flowminder Foundation, Forrest R. Stevens, Dana R. Thomson
##  
##  Permission is hereby granted, free of charge, to any person obtaining
##  a copy of this software and associated documentation files (the
##  "Software"), to deal in the Software without restriction, including
##  without limitation the rights to use, copy, modify, merge, publish,
##  distribute, sublicense, and/or sell copies of the Software, and to
##  permit persons to whom the Software is furnished to do so, subject to
##  the following conditions:
##  
##  The above copyright notice and this permission notice shall be
##  included in all copies or substantial portions of the Software.
##  
##  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
##  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
##  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
##  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
##  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
##  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
##  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



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
