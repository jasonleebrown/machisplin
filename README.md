![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/MACHISPLIN_LOGOv2.jpg?raw=true "Title")

# machisplin
An R package for interpolation of noisy multi-variate data through comprehensive statistical analyses using thin-plate-smoothing splines and machine learning ensembling.  This package is a free open-source machine learning analog to the expensive ANUSPLIN software. 

## To install open R and type the following:
```markdown
install.packages("devtools")
library(devtools)
install_github("jasonleebrown/machisplin")
```

This script interpolates noisy multi-variate data through machine learning ensembling using six algorithms: boosted regression trees (BRT), neural networks (NN); generalized additive model (GAM), multivariate adaptive regression splines (MARS), support vector machines (SVM) and random forests (RF). This function evaluates (via k-fold cross validation, where k=10) a method’s ability to predict the input data and ensembles of all combinations of the six algorithms weighting each from 0-1 and evaluting fit. The best model will have the lowest AICc (with the number of parameters in AICc calculation corresponding the number of models in ensemble). After the best model is determined, the function will run the ensemble on the full dataset. Then residuals will be calculated and interpolated using thin-plate-smoothing splines, which will secondarily correct the final ensemble model. This package is a free open-source machine learning analog to the expensive ANUSPLIN software. To output final R2 values, model weights, algorithm(s) used, and rasters for use in GIS; use the 'machisplin.write.geotiff' function. To output residuals use 'machisplin.write.residuals' and to output model loadings use 'machispline.write.loadings'.

Example 1
```markdown
######## EXAMPLE 1 ########
library(MACHISPLIN)

# Import a csv as shapefile:
Mydata<-read.delim("sampling.csv", sep=",", h=T)

#load rasters to use as high resolution co-variates for downscaling
ALT = raster("SRTM30m.tif")
SLOPE = raster("ln_slope.tif")
ASPECT = raster("aspect.tif")
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

Example 2
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
![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/Slide19.JPG?raw=true "Title")
![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/Slide20.JPG?raw=true "Title")
![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/Slide21.JPG?raw=true "Title")
![Alt text](https://raw.githubusercontent.com/jasonleebrown/machisplin/master/Slide26.JPG?raw=true "Title")


[![Analytics](https://ga-beacon.appspot.com/UA-136960917-1/machisplin)](https://github.com/igrigorik/ga-beacon)
[![Analytics](https://ga-beacon.appspot.com/UA-136933757-1/machuruku?pixel)](https://github.com/igrigorik/ga-beacon)
