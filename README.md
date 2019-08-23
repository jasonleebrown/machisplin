![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/MACHISPLIN_LOGOv2.jpg?raw=true "Title")

# machisplin
An R package for interpolation of noisy multi-variate data through comprehensive statistical analyses using thin-plate-smoothing splines and machine learning ensembling.  This package is a free open-source machine learning analog to the expensive ANUSPLIN software. 

## To install open R and type the following:
```markdown
install.packages("devtools")
library(devtools)
install_github("jasonleebrown/machisplin")
```

## Contact 
[We can help you sort out issues, maybe](https://www.jasonleebrown.org/get-in-touch)

## Input Data formats
To explore data format for input data, see:
```markdown
library(humboldt)
##data format that will be interpolated, each column is a different dataset
Mydata<-sampling

#rasters to use as high resolution co-variates for downscaling
ALT = raster(system.file("extdata", "alt.tif", package="MACHISPLIN"))
SLOPE = raster(system.file("extdata", "slope.tif", package="MACHISPLIN"))
TWI = raster(system.file("extdata", "TWI.tif", package="MACHISPLIN"))
```
If needed, see guide below to convert raster GIS data for use as in 'MACHISPLIN'


## Analyses in 'MACHISPLIN'
This script interpolates noisy multi-variate data through machine learning ensembling using six algorithms: boosted regression trees (BRT), neural networks (NN); generalized additive model (GAM), multivariate adaptive regression splines (MARS), support vector machines (SVM) and random forests (RF). This function evaluates (via k-fold cross validation, where k=10) a method’s ability to predict the input data and ensembles of all combinations of the six algorithms weighting each from 0-1 and evaluting fit.
The best model will have the lowest AICc (with the number of parameters in AICc calculation corresponding the number of models in ensemble). After the best model is determined, the function will run the ensemble on the full dataset. Then residuals will be calculated and interpolated using thin-plate-smoothing splines, which will secondarily correct the final ensemble model. This package is a free open-source machine learning analog to the expensive ANUSPLIN software. To output final R2 values, model weights, algorithm(s) used, and rasters for use in GIS; use the 'machisplin.write.geotiff' function. To output residuals use 'machisplin.write.residuals' and to output model loadings use 'machispline.write.loadings'.
![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/Slide20.JPG?raw=true "Title")
Overview of Process

![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/Slide21.JPG?raw=true "Title")
Details of Modeling using Hi-resolution Covariates 

![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/Slide26.JPG?raw=true "Title")
Example of the results (all with R2>0.99).

### Example 1 - uisng provide datasets
```markdown
library(MACHISPLIN)
library(raster)

##load spatial data with (coordinates named exactly as 'long' and 'lat') and any number of layers to downscale
data(sampling)
Mydata<-sampling

#load rasters to use as high resolution co-variates for downscaling
ALT = raster(system.file("extdata", "alt.tif", package="MACHISPLIN"))
SLOPE = raster(system.file("extdata", "slope.tif", package="MACHISPLIN"))
TWI = raster(system.file("extdata", "TWI.tif", package="MACHISPLIN"))

# function input: raster brick of covarites
raster_stack<-stack(ALT,SLOPE,TWI)

#run an ensemble machine learning thin plate spline 
interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2)

machisplin.write.geotiff(mltps.in=interp.rast)
machisplin.write.residuals(mltps.in=interp.rast)
machisplin.write.loadings(mltps.in=interp.rast)
```
### Example 2 - typical workflow
see below for help formating raster/environment data. 
```markdown
library(MACHISPLIN)

# Import a csv as shapefile:
Mydata<-read.delim("sampling.csv", sep=",", h=T)

#load rasters to use as high resolution co-variates for downscaling
ALT = raster("SRTM30m.tif")
SLOPE = raster("ln_slope.tif")
TWI = raster("TWI.tif")

# function input: raster brick of covarites
raster_stack<-stack(ALT,SLOPE,TWI,ASPECT)

#run an ensemble machine learning thin plate spline 
interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2, tps=FALSE)

#to plot results (change number to select different output raster)
plot(interp.rast[[1]]$final)

#to view residuals (change number to select different output raster)
interp.rast[[1]]$residuals

#to view model loadings
interp.rast[[1]]$var.imp

#to view other model parameters and other parameters
interp.rast[[1]]$summary
```

## Getting environmental data formatted for ‘MACHISPLIN’
You need two sets of datafiles.  1. the layers that will be interpolated and 2. higher resolution covariates that will be use to downscale interpolation layers 

### 1. Layers to be interpoloated 
This is a single data file where the first two columns are longitude and latitude (x and y) in that order.  The following columns represent the corresponding values of the data layers that will be interpoloated.   This can a single layer (=1 column) or a dozen (=12 columns).    

### 2. High resolution covariates
These need to be a series of high resolution raster combinded into a raster stack.   All rasters must be the same: resolution, projection and extent.  Typically these are microtopgraphic layers.  

### Importing rasters into R
```markdown
library(raster)

##rasters to downscale, here Tiff rasters are in base working directory
BIO1 = raster("bio1.tif")
BIO2 = raster("bio2.tif")
BIO12 = raster("bio12.tif")
BIO15= raster("bio15.tif")

#high resolution covariate rasters
ALT = raster("SRTM30m.tif")
SLOPE = raster("ln_slope.tif")
TWI = raster("TWI.tif")

```
### Clipping the extent of rasters in R
```markdown
##create raster stack to speed process up 
raster_interp_layers<-stack(BIO1, BIO2, BIO12, BIO15)
raster_covar_layers<-stack(ALT, SLOPE, TWI)

##define extent for which you will clip data to. These values will reflect the spatial area of analyses and the area where you will create high resolution variables.    
##I highly recommend determining these values from a GIS (ArcGIS or Google Earth).
##extent order input below is: xmin, xmax, ymin, yma
e.in <- extent(-160, 10, 30, 60)

##perform crop function to clip to input extent
interp_layers <- crop(raster_interp_layers, e.in)	
covar_layers <- crop(raster_covar_layers, e.in)	
```

### Convert input 'rasters to interpolated' to 'MACHISPLIN' format (continued from previous step)
```markdown
##convert one of the rasters to a point dataframe to sample.  Use any raster input.
env.points<-rasterToPoints(BIO1, fun=NULL, spatial=FALSE)

##subset only the x and y data
env.points<- env.points[,1:2]

##Extract values to points from rasters
RAST_VAL<-data.frame(extract(interp_layers, env.points))

##merge sampled data to input
InInterp<-cbind(env.sampling.res,RAST_VAL)

##save the file as '.csv' for future analyses 
write.csv(Env1, file = "InInterp_v1.csv")
```

Now use the 'InInterp' and 'raster_covar_layers' in the 'machisplin.mltps' function, for example:
'interp.rast<-machisplin.mltps(int.values=InInterp, covar.ras=raster_covar_layers, n.cores=2)'


# Need help with the high resolution topography data? 
1. Get a high resolution elevation model
-Download 30m or 90m elevation data and mosaic tile to single file in GIS
https://earthexplorer.usgs.gov/ (requires registration)

2. Use the high resolution elevation model to create topography layers:
Geomorphons:
 use Grass GIS (http://shorturl.at/cesZ1);
Topographic Wetness Index (TWI), slope and aspect:
 use SAGA GIS (http://www.saga-gis.org/)
 Note that you can also generate slope and aspect in ArcMap or R, but in these programs the creation of TWI is quite complex.  Thus its best just to use SAGA GIS in which all are done by simply clicking a button.


[![Analytics](https://ga-beacon.appspot.com/UA-136960917-1/machisplin)](https://github.com/igrigorik/ga-beacon)
[![Analytics](https://ga-beacon.appspot.com/UA-136933757-1/machuruku?pixel)](https://github.com/igrigorik/ga-beacon)
