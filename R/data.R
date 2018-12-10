#' Example Data: input table of values and rasters to downscale
#' 
#' @docType data
#' 
#' @usage data(to.interp)
#' 
#' @format An data frame with 813 rows and 4 variables (coordinates named exactly as 'long' and 'lat') and two layers to downscale:
#' \describe{
#'   \item{long}{longitude}#long,lat,bio_1,bio_12
#'   \item{lat}{latittude}
#'   \item{var1}{bio_1}
#'   \item{var2}{bio_12}
#' }
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' library(MACHISPLIN)
#' library(raster)
#' 
#' ##load spatial data with (coordinates named exactly as 'long' and 'lat') and any number of layers to downscale
#' data(sampling)
#' Mydata<-sampling
#' 
#' ## load high-resolution covariates rasters
#' 
#' ALT = raster(system.file("extdata", "alt.tif", package="MACHISPLIN"))
#' 
#' SLOPE = raster(system.file("extdata", "slope.tif", package="MACHISPLIN"))
#' 
#' TWI = raster(system.file("extdata", "TWI.tif", package="MACHISPLIN"))
#' 
#' ##function input: raster stack of covariates
#' raster_stack<-stack(ALT,SLOPE,TWI)

"example.dat"

