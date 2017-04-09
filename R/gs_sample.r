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



globalVariables(names=c("selected","id","grid_id","stratum","V1","sample_ids","psu_id","valid","pop","urban"))
#' @title GridSample sampling algorithm
#' @description
#' The \code{gs_sample} algorithm creates primary sampling units (PSUs) for multi-stage cluster household surveys based on gridded population data. Typical complex survey design is supported with input of a raster of population counts, a raster of urbanized areas, and a raster of study strata. Each of these rasters need to be in an identical projection and have an identical grid resolution. The algorithm first selects PSU seed cells with a probability proportionate to population size according to strata, rural-urban, and spatial parameters specified, then it optionally grows PSUs around the seed cells until a minimum population threshold is met in each PSU.
#' 
#' @details
#' A number of sampling features are optional. Oversampling in urban/rural areas, oversampling to be spatially representative, and stratification are not required. At a minimum, the user generates a simple random sample of PSUs in a study area by inputting a \code{population_raster}, defining the study area boundary as one stratum with \code{strata_raster}, defining the output shapefile parameters \code{output_path} and \code{sample_name}, and configuring the parameters \code{cfg_hh_per_stratum}, \code{cfg_hh_per_urban}, \code{cfg_hh_per_rural}, and \code{cfg_pop_per_psu}.  See the "Stratification", "Urban/rural domains", "Spatial sampling", and "PSU size and framework" sections for additional information. Note that all datasets are re-projected into WGS84 before the sampling process begins. A real-world example can be seen using the code \code{vignette("Rwanda")}, a vignette that replicates the sample design of the 2010 Rwanda DHS survey.
#'
#' @section Stratification:
#' To stratify the sample, define geographic strata boundaries with \code{strata_raster}, and specify the sample size per strata with \code{cfg_hh_per_stratum}. For example, if a national survey will sample 10,000 households from 5 provinces, then \code{cfg_hh_per_stratum = 2000}. The parameter \code{cfg_hh_per_stratum} is the minimum sample size to generate representative population statistics. In some surveys, strata follow urban/rural boundaries within administrative units. If this is the case, then \code{strata_raster} should include the boundaries of urban and rural sampling areas within each administrative area, and \code{cfg_hh_per_stratum} should reflect the correct sample size per stratum - for example, a national sample of 10,000 households from each urban and rural areas in 5 provinces would have \code{cfg_hh_per_stratum = 1000}.  
#' @section Urban/rural domains:
#' If urban/rural populations are not part of the stratification scheme, then they are often treated as a sub-domain. Sub-domains represent important sub-populations for which representative statistics are generated from the survey data, and thus each sub-domain (at the national-level) should meet the minimum sample size specified for each stratum. If either the urban/rural sub-domain does not include enough households to generate population statistics with the desired precision, then extra PSUs are oversampled in the smaller sub-domain. To implement this step with \code{gs_sample}, set \code{cfg_sample_rururb = 1}. In practice, rural areas are often more difficult and expensive to visit, and thus a greater number of households might be sampled from rural PSUs than urban PSUs. This is why the user may specify different numbers of households to be sampled from each urban PSUs (\code{cfg_hh_per_urban}) and rural PSUs (\code{cfg_hh_per_rural}); if the same number of households will be sampled from all PSUs, then configure both of these parameters with the same value. Note, the number of PSUs that will be generated in each stratum is \code{cfg_hh_per_stratum} divided by some number between \code{cfg_hh_per_urban} and \code{cfg_hh_per_rural}.
#' @section Spatial sampling:
#' To select a sample that is both representative of the population and of space, set \code{cfg_sample_spatial = 1} and specify \code{cfg_sample_spatial_scale}, the spatial scale at which the sample should be representative. The spatial scale should be meaningful; for example, it will facilitate small area estimates with limited statistical error for administrative units that are smaller than the stratification units. Determining an appropriate spatial scale might take trial and error. If the study area has large regions of sparse population, a typical non-spatially representative sample will follow the population distribution and have large areas without a PSU. In this case, the user might need to increase the spatial resolution \code{cfg_sample_spatial_scale} of the sample, or force the algorithm to generate more PSUs in each stratum by increasing \code{cfg_hh_per_stratum} and/or reducing \code{cfg_hh_per_urban} and \code{cfg_hh_per_rural}. 
#' @section PSU size and fieldwork:
#' Four additional parameters can be configured to deal with idiosyncrasies of gridded population data and improve feasibility of fieldwork. The user can set a maximum geographic size of PSU in kilometres squared, \code{cfg_max_psu_size}. We recommend choosing a size that can feasibly be visited by a field team on foot during one day. The user might also specify which cells are included in the sample frame with \code{cfg_min_pop_per_cell}. Selection of a sensible value is highly dependent on the gridded population dataset being used, and the scale of the input data (e.g. 200m X 200m grid cells). The cell size of the output raster can be specified with \code{cfg_desired_cell_size}. Gridded population datasets generated from old population figures or old covariates may be inaccurate at a very local scale (e.g. 100m X 100m cells), but will generally increase in accuracy as cells are aggregated (e.g. 300m X 300m cells). Finally, the PSU growth portion of the algorithm can be switched off by setting \code{cfg_psu_growth = FALSE} resulting in a sample of single grid cells (and their centroids).

#' @import raster rgdal geosphere rgeos maptools spatstat.utils
#' @importFrom utils flush.console 
#' @importFrom data.table data.table setkey
#' @importFrom spatstat dirichlet
#' @importFrom methods as
#' @param population_raster Raster* layer. Input gridded population dataset to use as sample frame. Values should be number of people in each pixel as a whole number or decimal value.
#' @param strata_raster Raster* layer. Raster that defines the stratum numeric ID of each pixel. Generally created by rasterizing a shapefile of polygons that define strata.
#' @param urban_raster Raster* layer. Raster of urbanized areas where a cell value of 1 indicates urban cells and 0 indicates rural cells.
#' @param cfg_hh_per_stratum numeric. Target household sample size per stratum. In a non-stratified sample, this is the total sample size of households. In a stratified sample, this is the household sample size per stratum.
#' @param cfg_hh_per_urban numeric. Number of households expected to be selected per urban PSU during survey fieldwork.
#' @param cfg_hh_per_rural numeric. Number of households expected to be selected per rural PSU during survey fieldwork.
#' @param cfg_pop_per_psu numeric. Minimum population per PSU (e.g. 500 persons).
#' @param cfg_sample_rururb logical. A flag to oversample rural/urban areas if one domain does not meet the target sample size per stratum. Default is \code{FALSE}.
#' @param cfg_sample_spatial logical. A flag to oversample in space ensuring that at least one PSU is selected within each "coarse grid" cell with cell size defined by the user. Default is \code{FALSE}.
#' @param cfg_sample_spatial_scale If \code{cfg_sample_spatial == TRUE}, this defines the size in kilometres squared (e.g. 20 for 20km X 20km) of each course grid cell where the algorithm will ensure at least one PSU is located in each coarse grid cell.
#' @param cfg_desired_cell_size numeric. Desired cell size in kilometres squared (e.g. 0.4 for 400m X 400m) for output raster of PSUs. Defaults to NA, which yields an output raster at the same resolution as population_raster.
#' @param cfg_max_psu_size numeric. Maximum allowed geographic size of a given PSU in kilometres squared (e.g. 5 for PSUs smaller than 5km X 5km). Defaults to infinity.
#' @param cfg_min_pop_per_cell numeric. Minimum population in a raster cell required for it to be considered for sampling. Cells with less than this value will be excluded from the sample. Defaults to 0, therefore including all cells.
#' @param cfg_psu_growth logical. Determines whether to grow PSUs until either there are no available cells or each PSU covers a population defined by \code{cfg_pop_per_psu}.
#' @param output_path character. Output path and folder name.
#' @param sample_name character. Name of output PSU shapefile.
#' 
#' @return Shapefile of household survey primary sampling unit (PSU) boundaries
#' @examples 
#' require(raster)
#' poprast <- raster(ncols=100,nrows=100,xmx=10,xmn=9,ymn=9,ymx=10,
#' crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
#' vals=runif(10000,0,100))
#' stratarast <- raster(ncols=100,nrows=100,xmx=10,xmn=9,ymn=9,ymx=10,
#' crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
#' vals=c(rep(1,times=5000),rep(2,times=5000)))
#' urbanrast <- poprast > 25
#' example_1 <- gs_sample(population_raster = poprast, 
#'          strata_raster = stratarast, 
#'          urban_raster = urbanrast, 
#'          cfg_hh_per_stratum = 800,
#'          cfg_hh_per_urban = 20, 
#'          cfg_hh_per_rural = 20, 
#'          cfg_pop_per_psu = 500,
#'          cfg_sample_rururb = TRUE, 
#'          cfg_sample_spatial = FALSE, 
#'          cfg_sample_spatial_scale = 100,
#'          cfg_desired_cell_size = NA,
#'          cfg_max_psu_size = 5,
#'          cfg_min_pop_per_cell = 0.01,
#'          output_path=tempdir(),
#'          sample_name="Example")
#' plot(example_1)
#' 
#' 
#' #### Example two is the identical, except PSUs aren't grown, 
#' #### so the shapefile returned includes a single grid cell for each PSU.
#' 
#' example_2 <- gs_sample(population_raster = poprast, 
#'          strata_raster = stratarast, 
#'          urban_raster = urbanrast, 
#'          cfg_hh_per_stratum = 800,
#'          cfg_hh_per_urban = 20, 
#'          cfg_hh_per_rural = 20, 
#'          cfg_pop_per_psu = 500,
#'          cfg_sample_rururb = TRUE, 
#'          cfg_sample_spatial = FALSE, 
#'          cfg_sample_spatial_scale = 100,
#'          cfg_desired_cell_size = NA,
#'          cfg_max_psu_size = 5,
#'          cfg_min_pop_per_cell = 0.01,
#'          cfg_psu_growth = FALSE,
#'          output_path=tempdir(),
#'          sample_name="Example_without_growth")
#' plot(example_2)
#' @export
#' 
gs_sample <- function(population_raster, 
                      strata_raster, 
                      urban_raster, 
                      cfg_hh_per_stratum,        
                      cfg_hh_per_urban, 
                      cfg_hh_per_rural, 
                      cfg_pop_per_psu,           
                      cfg_sample_rururb = FALSE, 
                      cfg_sample_spatial = FALSE, 
                      cfg_sample_spatial_scale = NA,
		                  cfg_desired_cell_size =NA,
                      cfg_max_psu_size = Inf,
                      cfg_min_pop_per_cell = 0,
		                  cfg_psu_growth = TRUE,
		                  output_path,
		                  sample_name){
  ## Project rasters in WGS84 if not
  if (sp::proj4string(population_raster) != "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"){
  population_raster <- projectRaster(population_raster, crs=sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
  if (sp::proj4string(urban_raster) != "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"){
    urban_raster <- projectRaster(urban_raster, crs=sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
  
  if (sp::proj4string(strata_raster) != "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"){
    strata_raster <- projectRaster(strata_raster, crs=sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
  
  ## Aggregate to meet desired resolution if lower than current
   if (!is.na(cfg_desired_cell_size)){
    
    midcell <- c(xFromCol(population_raster,col=floor(ncol(population_raster)/2)),yFromRow(population_raster,row=nrow(population_raster)/2))
    nextcell <- c(xFromCol(population_raster,col=floor(ncol(population_raster)/2)+1),yFromRow(population_raster,row=nrow(population_raster)/2))
    
    cellsize <- distGeo(midcell,nextcell)/1000
    if (cellsize < cfg_desired_cell_size){
      population_raster <- aggregate(population_raster,fact=round(sqrt(cfg_desired_cell_size)/cellsize),fun=sum,na.rm=T)
      strata_raster <- resample(strata_raster,population_raster,method="ngb")
      urban_raster <- resample(urban_raster,population_raster,method="ngb")
      
    }
    if (cellsize >= cfg_desired_cell_size){
      print("Resolution already as high or higher than desired.")
    }
  }
  midcell <- c(xFromCol(population_raster, col = floor(ncol(population_raster)/2)),yFromRow(population_raster, row = nrow(population_raster)/2))
  nextcellx <- c(xFromCol(population_raster, col = floor(ncol(population_raster)/2)+1), yFromRow(population_raster, row = nrow(population_raster)/2))
  nextcelly <- c(xFromCol(population_raster,col = floor(ncol(population_raster)/2)), yFromRow(population_raster, row = nrow(population_raster)/2+1))
  
  cellsizex <- distGeo(midcell,nextcellx)/1000 #convert to kilometers
  cellsizey <- distGeo(midcell,nextcelly)/1000
  
  cellkmsq <- cellsizex * cellsizey
  if (is.na(cfg_sample_spatial_scale)==FALSE){
  spatialfac=(cfg_sample_spatial_scale)/mean(c(cellsizex,cellsizey))
  }
  ## Convert maximum PSU size (in square kilometers) into a maximum possible number of cells
  maxcells <- cfg_max_psu_size / cellkmsq
	print("Beginning gs_sample() Process:")
	flush.console()
	
	
	##	Compile vectors of values:
	urban_vec <- urban_raster[]
	pop_vec <- population_raster[]
	strata_vec <- strata_raster[]

	##	Create an index of cells, and vectors with no NAs from our data:
	raster_index <- 1:ncell(strata_raster)
	raster_index_valid_boo <- (!is.na(strata_vec) & !is.na(urban_vec) & !is.na(pop_vec) & pop_vec > cfg_min_pop_per_cell) 
	raster_index_valid <- raster_index[ raster_index_valid_boo ]

	##	Create sampled_index that will be values for our output, 0 for 
	##		unsampled, 1 for sampled:
	sampled <- 0 * raster_index_valid_boo

	##	Compile unique lists of valid Strata IDs:
	strata_ids <- unique( strata_vec[raster_index_valid] ) 
	print("Randomly Stratifying PSU Cell IDs:")
	flush.console()


    for (stratum_id in strata_ids) {
      print(paste("Processing stratum ID", stratum_id, "..."))
      flush.console()
      ##	Create index for rasters for current strata, with valid values:
      index <- raster_index[ raster_index_valid_boo & (strata_vec == stratum_id) ]
      initsampledpsu <- cfg_hh_per_stratum / min(c(cfg_hh_per_urban,cfg_hh_per_rural),na.rm=T)
      psu_data <- data.table("id"=index, "pop"=pop_vec[index],"urban"=urban_vec[index],"selected"=0)
      ## Calculate number of cells to skip for serpentine sampling
      psu_interval <- floor(sum(psu_data$pop) /  initsampledpsu)
      
      psu_data$pop_cum <- cumsum(psu_data$pop)
      psu_data$urbanpop <- (psu_data$urban == 0) * cfg_hh_per_rural + (psu_data$urban == 1) * cfg_hh_per_urban
      psu_data$ruralpop <- (psu_data$urban == 1) * cfg_hh_per_rural + (psu_data$urban == 0) * cfg_hh_per_urban
      
      ## Choose first cell to sample
      initpop=sample(psu_interval,1)
      psu_data$selected[which.min(abs(psu_data$pop_cum - initpop))] <- 1
      
      ## Skip based on psu_interval then select cells
      nextpop <- initpop + psu_interval
      while (nextpop<max(psu_data$pop_cum)){
        psu_data$selected[which.min(abs(psu_data$pop_cum - nextpop))] <- 1
        nextpop <- nextpop + psu_interval
      }
      ## Calculate maximum number of households sampled given the cells that have been chosen
      psu_data$maxsample <- max(psu_data$urbanpop,psu_data$ruralpop)
      total_stratum <- sum(psu_data[selected==1]$maxsample)
      
      ## Remove and add cells as necessary based on how many will be sampled per stratum, compared to the desired cfg_hh_per_stratum
      if (total_stratum > cfg_hh_per_stratum){
        popremove <- floor(total_stratum - cfg_hh_per_stratum)
        psu_seeds <- psu_data[selected == 1]
        psu_seeds <- psu_seeds[sample(1:dim(psu_seeds)[1]), ]
        psu_seeds$popcum <- cumsum(psu_seeds$maxsample)
        psu_data$selected[is.element(psu_data$id, psu_seeds$id[which(psu_seeds$urbanpopcum < popremove)])] <- 0
                }
      if (total_stratum < cfg_hh_per_stratum){
        while (total_stratum < cfg_hh_per_stratum){
          newid <- sample(psu_data[selected == 0]$id, 1)
          total_stratum <- total_stratum + (urban_vec[newid] == 0) * cfg_hh_per_rural + (urban_vec[newid] == 1) * cfg_hh_per_urban
          psu_data[id == newid]$selected <- 1
        }
      }
      
    ## set all selected cells to 1 in the "sampled" vector
      sampled[ psu_data[selected==1]$id ] <- 1
    }
	##	Next we need to revisit where we have 1 in our sampled values, and 
	##		confirm	that we match our rural/urban ratio if 
	##		cfg_sample_rururb == TRUE:

	if (cfg_sample_rururb == TRUE) {
		print("Checking urban and rural domains meet stratum HH sample size")
		flush.console()

		##	Extract the rural and urban values that are currently in our sample:
		urban_vec_sampled <- urban_vec[sampled == 1]
	
		rural_count <- length(which(urban_vec_sampled == 0))
		urban_count <- length(which(urban_vec_sampled == 1))
	
		##	Double check to see if we are missing PSUs of a certain type:
		reassign <- NA
		choice <- NA
    num_switch <- 0
		##	Find out which class we need to increase, and which we need to 
		##		decrease:
		if (rural_count * cfg_hh_per_rural < cfg_hh_per_stratum) {
			##	We need more rural in our sample:
			choice <- 1
			reassign <- 0
			
			##	Number to reassign:
			num_switch <- round(cfg_hh_per_stratum/cfg_hh_per_rural) - rural_count
		}
		if (urban_count * cfg_hh_per_urban < cfg_hh_per_stratum) {
			##	We need more urban in our sample:
			choice <- 0
			reassign <- 1
			
			##	Number to reassign:
			num_switch <- round(cfg_hh_per_stratum/cfg_hh_per_urban) - urban_count
		}
		if (rural_count * cfg_hh_per_rural < cfg_hh_per_stratum & urban_count*cfg_hh_per_urban < cfg_hh_per_stratum) {
		print("Both are below threshold! Cannot reassign to meet desired sample size")
    reassign <- NA
    choice <- NA
		  		}
		if (is.na(reassign)==TRUE){print("Not reassigning urban/rural")}
    
    if (is.na(reassign)==FALSE){
      
		reassign_index <- raster_index[ (sampled == 0) & (urban_vec == reassign) & raster_index_valid_boo ]
		choice_index <- raster_index[ (sampled == 1) & (urban_vec == choice) & raster_index_valid_boo ]
		if (length(reassign_index)==0 | length(choice_index)==0){
		  print("Unable to reassign urban/rural; no available cells")
		  
		  
		}
		## Choose cells to reassign
		if (length(reassign_index) > 0 & length(choice_index) > 0){
		reassign_index = reassign_index[is.element(strata_vec[reassign_index], unique(strata_vec[choice_index]))]
		reassign_ids <- sample(reassign_index, num_switch)
		
		reassign_strata <- strata_vec[reassign_ids]
		strata_numbers <- table(reassign_strata)
		strata_numbers_ids <- as.numeric(names(strata_numbers))

		if (length(reassign_strata) > 0) {
			for (i in 1:length(strata_numbers_ids)) {
				stratum_id <- strata_numbers_ids[i]
				stratum_number <- strata_numbers[[i]] 
				if (choice == 1) {
					stratum_number <- floor((cfg_hh_per_urban/cfg_hh_per_rural)*stratum_number)
				} else {
					stratum_number <- floor((cfg_hh_per_rural/cfg_hh_per_urban)*stratum_number)
				}
				
				print(paste("Reassignment for stratum", strata_ids[i], "of", stratum_number, "cells processed"))
				choice_index <- raster_index[ (sampled == 1) & (urban_vec == choice) & strata_vec == stratum_id ]
				##	If the number of samples to convert is more than 90% of the 
				##		total number then only convert up to that amount:
				if (stratum_number > round(0.90 * length(choice_index))) {
					stratum_number <- round(0.90 * length(choice_index))
				}

				choice_ids <- sample(choice_index, stratum_number)

				sampled[ choice_ids ] <- 0
			}

			##	Finally, convert our reassigned cells to 1 in the sampled vector:
			sampled[ reassign_ids ] <- 1
		}
  }
    }
	}
  
  

	if (cfg_sample_spatial == TRUE) {
		print("Stratifying by Spatial Grid:")
		flush.console()

		##	First create a raster with cell indices from our valid cell 
		##		locations:
		valid_cells <- raster_index
		valid_cells[-1*raster_index_valid] <- 0
		valid_raster <- strata_raster
		valid_raster <- setValues(valid_raster, valid_cells)
		## Change the aggregated cell parameter from square kilometers to numbers of pixels
		valid_raster <- aggregate(valid_raster, fac = spatialfac, fun = max)
		valid_raster <- resample(valid_raster, strata_raster, method = "ngb")
		
		valid_vec <- valid_raster[]
    

		##	Compile unique lists of valid aggregated grid IDs:
		grid_ids <- unique( valid_vec[raster_index_valid] )

		##	Compile counts of sampled points by grid ids:
		raster_data <- data.table("grid_id" = valid_vec, "sampled" = as.numeric(sampled), "stratum" = strata_vec, "raster_index" = raster_index)
		data.table::setkey(raster_data, grid_id, stratum, sampled)
		grid_table <- raster_data[, sum(sampled, na.rm=TRUE), by=grid_id]

		##	Compile majority value of the stratum id for each cell:
		strata_table <- raster_data[, gs_mode(stratum), by=grid_id]

		##	For each empty cell, randomly select a new cell within it to sample
		##		and remove one randomly from the stratum that makes up the majority
		##		of the cell:
		for (i in 1:length(strata_table[, grid_id])) {
			##	Check to see if the grid is empty and the dominant stratum is not 
			##		zero:
			if (grid_table[, V1][i] == 0 & !is.na(strata_table[, V1][i])) {

				print(paste("Reassignment for grid", grid_table[, grid_id][i], "and stratum", strata_table[, V1][i]))
				flush.console()


				##	Select unsampled raster indices for this grid cell:
				cell_raster_indices <- raster_data[list( strata_table[, grid_id][i], strata_table[, V1][i], 0)]
			
				##	Select sampled raster indices for this stratum:
				stratum_raster_indices <- raster_data[stratum == strata_table[, V1][i] & sampled == 1]
				
				##	Randomly choose a valid cell from this grid location to sample 
				##		and one from the sampled locations in the majority stratum to 
				##		unsample:
				compiled_cell_indices <- cell_raster_indices[, raster_index]
				compiled_cell_indices <- compiled_cell_indices[ compiled_cell_indices %in% raster_index[ (sampled == 0) & raster_index_valid_boo ]]
				
				if (length(compiled_cell_indices) > 0 & length(stratum_raster_indices[,raster_index]) > 0) {
					cell_raster_index <- sample(compiled_cell_indices, 1)

					##	Make sure we only sample from indices matching the rural/urban
					##		designation we just chose randomly:
					stratum_raster_indices <- raster_data[stratum == strata_table[,V1][i] & sampled == 1 & urban_vec == urban_vec[cell_raster_index]]
					stratum_raster_index <- sample(stratum_raster_indices[,raster_index], 1)

					sampled[ cell_raster_index ] <- 1
					sampled[ stratum_raster_index ] <- 0
				} else {
					print("No valid cells for this stratum left for assignment or reassignment!  Skipping...")
					flush.console()
				}
			}
		}
	}
  ## Convert sampled vector into a raster
	sampled_image <- strata_raster
	sampled_image <- setValues(sampled_image, sampled)
	
	##	Second Stage Sampling:

	print("Begin Second Stage Sampling, Growing PSUs if requested:")
	flush.console()

	##	Assign an ID to each sampled cell randomly:
	sampled_ids <- sample( 1:sum(sampled) )
	sampled_image[ sampled_image[] != 0 ] <- sampled_ids
	sampled_image[ sampled_image[] == 0 ] <- NA
	cellids=which(is.na(sampled_image[])==F)
	xydat=data.frame(sampleid = sampled_image[cellids],x=xFromCell(sampled_image,cell=cellids),y=yFromCell(sampled_image,cell=cellids))
	
	##	Create data.table of our data for sorting and quick accesss:
	raster_data <- data.table("sampled" = sampled, "sample_ids" = sampled_image[], "stratum" = strata_vec, "raster_index" = raster_index)
	setkey(raster_data, sample_ids)

	index_raster <- sampled_image
	index_raster[] <- raster_index

	##	Create vectors of x/y coordinates from our raster object for use
	##		in creating a distance from focal PSU point within the loop:
	xCoords <- xFromCell(index_raster, index_raster[])
	yCoords <- yFromCell(index_raster, index_raster[])

	##	Initialize a new image based on our sampled image for storing our
	##		expanded PSU areas:
	psu_image <- sampled_image
	psu_image[ is.na(psu_image[]) ] <- 0

  if (cfg_psu_growth == TRUE){


	psu_ids <- sampled_ids
	psu_sampled_ids <- sampled_ids
	psu_voronoi_done <- sampled_ids
	
	## Calculate statistics for each PSU (number of households, number of cells, etc, and a vector indicating if the Voronoi polygons are filled)
	psu_sampled_sum <- population_raster[][raster_data[list(psu_sampled_ids),raster_index]]
	psu_sampled_cell_count <- psu_sampled_sum*0 + 1
	psu_voronoi_done <- sampled_ids
	
	psu_index <- sp::SpatialPointsDataFrame(coords=xyFromCell(index_raster, raster_data[sampled==1]$raster_index),data=as.data.frame(raster_data[sampled==1]),proj4string=sp::CRS(sp::proj4string(population_raster)))
	## Create Voronoi polygons around each seed PSU cell
	psuvoronoi=as(dirichlet(spatstat::as.ppp(psu_index@coords,W=extent(population_raster)[1:4])),"SpatialPolygons")
	psuvoronoi$sample_ids=psu_index$sample_ids
	psuvoronoi$stratum=psu_index$stratum
	voronoirast=rasterize(psuvoronoi,population_raster,field="sample_ids")
	passes <- 0
	while (sum(psu_sampled_ids != 0) > 0) {

		print(paste("PSUs to process:",sum(psu_sampled_ids != 0)))
		flush.console()

		##	Randomly shuffle the remaining IDs so we don't process them in 
		##		the same order:
		shuffle_index <- sample(1:length(psu_sampled_ids))
  
		psu_ids <- psu_ids[shuffle_index]
		psu_sampled_ids <- psu_sampled_ids[shuffle_index]
		psu_sampled_sum <- psu_sampled_sum[shuffle_index]
		psu_sampled_cell_count <- psu_sampled_cell_count[shuffle_index]
		psu_voronoi_done <- psu_voronoi_done[shuffle_index]
		
		##	Print our current active PSUs:
		print("Active PSUs:")
		print(psu_sampled_ids[])
		print("Active PSU Population Sums:")
		print(psu_sampled_sum[])
		print("Active PSU Cell Count:")
		print(psu_sampled_cell_count)
		
		##	Remove any PSUs that already exceed our population sum threshold:
		psu_sampled_ids[psu_sampled_sum[] > cfg_pop_per_psu] <- 0
		psu_sampled_ids[psu_sampled_cell_count > maxcells] <- 0
		psu_voronoi_done[psu_sampled_sum[] > cfg_pop_per_psu] <- 0
		psu_voronoi_done[psu_sampled_cell_count > maxcells] <- 0
		
	  if (length( which( psu_voronoi_done > 0) > 0)){
		

		for (sample_id in psu_voronoi_done[psu_voronoi_done > 0]) {
		  
			print(paste("Processing",sample_id,"..."))
			flush.console()
    
			##	Calculate distances from the current cell for entire image:
			x_index <- mean(xFromCell(index_raster,which(psu_image[] == sample_id)))
			y_index <- mean(yFromCell(index_raster,which(psu_image[] == sample_id)))
			
			cell_dist <- ((x_index - xCoords)^2 + (y_index - yCoords)^2)^0.5

			##	Create an index that only pulls values for those cells within our
			##		current stratum, and for which we have not already selected for
			##		use in a PSU:
			stratum_id <- raster_data[list(sample_id),stratum]
			stratum_indices <- raster_data[stratum==stratum_id,raster_index]
			indices <- stratum_indices[ psu_image[stratum_indices] == 0 & raster_index_valid_boo[stratum_indices] == T & voronoirast[stratum_indices] == sample_id]

			if (length(indices) == 0) {
				##	If there are no more valid indices left in the strata for 
				##		assignment then we drop the ID from the psu_sampled_ids list:
				psu_sampled_ids[psu_sampled_ids == sample_id] <- 0
			  print("No Voronoi-bounded cells left!")
			} else {
				
				ordered_data <- data.frame(indices, "pop"=pop_vec[indices], "dist"=cell_dist[indices])
				ordered_data <- ordered_data[order(ordered_data$dist, -(ordered_data$pop),decreasing=F),]

				ordered_data$pop_cum <- cumsum(ordered_data$pop) + psu_sampled_sum[psu_sampled_ids == sample_id]

				selection <- ordered_data$pop_cum <= cfg_pop_per_psu 
				selection[which(selection==FALSE)[1]] <- TRUE
				
				indices_to_add <- ordered_data[ selection, "indices" ]
				max_index <- length(indices_to_add)

				
				
				if (length(indices_to_add) == 0) {
					max_index <- 1
					indices_to_add[max_index] <- ordered_data[max_index,1]
				}
					
				if (sum(is.na(indices_to_add) > 0)) {
					indices_to_add <- indices_to_add[!is.na(indices_to_add)]
				}


				##	Update our psu_image based on new sample indices:
				psu_image[ indices_to_add ]  <- sample_id


				##	Adjust the total sample_id to match the total population:
				sum_population_added <- sum(pop_vec[indices_to_add])
				sum_population_for_sample <- sum(pop_vec[raster_index[psu_image[]==sample_id]])
				sum_cells_added <- length(indices_to_add)
				sum_cells <- sum(psu_image[] == sample_id)
				
				psu_sampled_sum[psu_sampled_ids == sample_id] <- sum_population_for_sample
				psu_sampled_cell_count[psu_sampled_ids == sample_id] <- sum_cells
				
				
				print(paste("... cells added:", sum_cells_added))
				print(paste("... population added:", sum_population_added))
				print(paste("... population total:", sum_population_for_sample))
				

				##	Check to see if we have reached our maximum configured 
				##		population or if we could not add any more cells, then 
				##		drop the current sample_id from the list to process:
				if ((sum_cells_added == 0) | (sum_population_for_sample > cfg_pop_per_psu)) {
					psu_sampled_ids[psu_sampled_ids == sample_id] <- 0
					psu_voronoi_done[psu_voronoi_done==sample_id]<-0
				}
			}
		}
		passes <- passes + 1
		}
		}
  }

	##	Last let's convert our 0's to NA and save the image off:
	psu_image[ psu_image[] == 0 ] <- NA

	print("Converting PSU Cells to Polygons and Summarizing Metrics:")

	##	Convert our PSUs to polygons:
	psu_polygons <- rasterToPolygons(psu_image,na.rm=T,dissolve=T)
	names(psu_polygons) <- "psu_id"

	#	Create data.table of our data for sorting and quick accesss of
	#		per-stratum information:
	raster_data <- data.table("valid" = raster_index_valid_boo, "sampled" = sampled, "sample_ids" = sampled_image[], "psu_id" = psu_image[], "stratum" = strata_vec, "pop" = as.numeric(population_raster[]), "urban" = urban_vec, "raster_index" = raster_index)

	setkey(raster_data, psu_id)

	#	Calculate, total, rural (urban==0) and urban (urban==1) populations
	#		per stratum, as well as total number of cells:
	stratum_pop <- raster_data[valid == TRUE, list("stratum"=max(stratum), "stratum_pop"=sum(pop, na.rm=TRUE), "stratum_rural_pop"=sum((urban == 0)*pop, na.rm=TRUE), "stratum_urban_pop"=sum((urban == 1)*pop, na.rm=TRUE), "stratum_cells"=length(raster_index)), by=stratum]
	setkey(stratum_pop, stratum)

	psu_pop <- raster_data[valid == TRUE & !is.na(psu_id), list("stratum"=max(stratum), "psu_pop"=sum(pop, na.rm=TRUE), "psu_rural_pop"=sum((urban == 0)*pop, na.rm=TRUE), "psu_urban_pop"=sum((urban == 1)*pop, na.rm=TRUE), "psu_cells"=length(raster_index)), by=psu_id]
	setkey(psu_pop, stratum)

	merge_pop <- psu_pop[list(stratum_pop)]
	setkey(merge_pop, psu_id)
	
	fields <- c("stratum", "psu_id", "psu_pop", "psu_rural_pop", "psu_urban_pop", "psu_cells", "stratum_pop", "stratum_rural_pop", "stratum_urban_pop", "stratum_cells")

	merge_pop <- as.data.frame(merge_pop)[,-7]
	psu_polygons@data = merge_pop[match(psu_polygons@data[,"psu_id"], merge_pop[,"psu_id"]),]
	psu_polygons = merge(psu_polygons,xydat,by.x="psu_id",by.y="sampleid")
  
	names(psu_polygons) <- c("psu_id", "stratum", "psu_pop", "psu_r_pop", "psu_u_pop", "psu_cells", "str_pop", "str_r_pop", "str_u_pop", "str_cells","xCent","yCent")

	writeOGR(psu_polygons, dsn=output_path, layer=sample_name, driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)

	return(psu_polygons)
}
