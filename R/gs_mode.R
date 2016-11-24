#' @author Forrest R. Stevens, \email{forrest.stevens@louisville.edu}
#' @title Most common stratum calculator
#' 
#' @description
#' For each cell in the coarse user-defined grid (only specified if \code{cfg_sample_spatial == TRUE}), this function calculates the stratum that occurs most often within each cell.
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
