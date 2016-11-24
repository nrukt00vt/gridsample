#' @author Forrest R. Stevens, \email{forrest.stevens@louisville.edu}
#' @title Rasterize polygon layer
#' 
#' @description
#' This function creates a raster layer that adopts values from a defined field in a polygon layer, using \code{rasterize} from the raster package.  This function also converts values to binary if desired, where all zero values are recorded as zero, and all non-zero values are recorded as one.  This function also saves the output raster in the working directory.
#'
#' 
#' @param input_features SpatialPolygons* object. Name of input shapefile layer. Should be a SpatialPolygons object.
#' @param output_raster Character. Desired name of output raster layer.
#' @param template_raster Raster* object. Raster with desired characteristics (resolution, extent) of output raster.
#' @param binary logical. If \code{TRUE}, any non-zero values will be converted to one.
#' @param field character. Name of variable that output raster should inherit.
#' @param overwrite logical. Defines whether to overwrite if \code{output_raster} already exists.
#' @param format character. Desired format of output raster file.

#' @return Vector of values representing the stratum that occurs most often within a given subset of the raster.
#' 
#' 
#' @export
gs_rasterize <- function(input_features, output_raster, template_raster, binary=FALSE, field="ID", overwrite=FALSE, format="GTiff") {
	if (file.exists(output_raster) & overwrite == FALSE) {
		strata_raster <- raster(output_raster)
	} else {
		strata_raster <- rasterize(input_features, template_raster, field=field)
		if (binary == TRUE) {
			strata_raster <- is.na(strata_raster)*0 + !is.na(strata_raster)*1
		}
		strata_raster <- writeRaster(strata_raster, filename=output_raster, format=format)
	}
	return(strata_raster)
}
