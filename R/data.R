#' Example Data: input table of values and rasters to downscale
#' 
#' @docType data
#' 
#' @usage data(sampling)
#' 
#' @format An data frame with 813 rows and 4 variables (oordinates named exactly as 'long' and 'lat') and two layers to downscale:
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
#' library(machisplin)
#' library(raster)
#' 
#' ##load spatial data with (coordinates named exactly as 'long' and 'lat') and any number of layers to downscale
#' data(sampling)
#' 
#' ## load environmental variables for all sites of the study area 2 (env2). Column names should be x,y,X1,X2,...,Xn)
#' 
#' ALT = raster(system.file("extdata", "alt.tif", package="machisplin"))
#' 
#' SLOPE = raster(system.file("extdata", "slope.tif", package="machisplin"))
#' 
#' TWI = raster(system.file("extdata", "TWI.tif", package="machisplin"))
#' 
#' ##function input: raster stack of covariates
#' raster_stack<-stack(ALT,SLOPE,TWI)

"sampling"

