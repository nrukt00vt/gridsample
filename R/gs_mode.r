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



#' @author Forrest R. Stevens, \email{forrest.stevens@louisville.edu}
#' @title Most common stratum calculator
#' 
#' @description
#' For each cell in the user-defined "coarse grid" used to select spatially-representative samples (i.e. when \code{cfg_sample_spatial == TRUE}), this function calculates the stratum that occurs most often within each course grid cell.
#'
#' 
#' @param rast \code{data.table} object. This \code{data.table} where each cell that lies within a larger grid cell is represented as a row. For each row, the variable  \code{grid_id} is the ID of the cell from the coarser grid,  \code{sampled} denotes whether a cell has been sampled,  \code{stratum} defines the stratum each cell lies within, and  \code{raster_index} is a unique value for each cell in the raster.
#' 
#' @return Vector of values representing the stratum that occurs most often within a given subset of the raster.
#' 
#' 
#' @export
gs_mode <- function(rast) {
	ux <- unique(rast)

	return(ux[which.max(tabulate(match(rast, ux)))])
}
