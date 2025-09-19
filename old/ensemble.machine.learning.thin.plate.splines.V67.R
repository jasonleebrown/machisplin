#'
#import nlme
#snowfall, snow, 

#add features to deal with empty grid for tilted TPS = grid= zero
##################################################################################################
##################################################################################################
#' Machine Learning Ensemble & Thin-Plate-Spline Interpolation
#' @param int.values An data frame with the first two columns as coordinates of interpolated site named exactly as 'long' and 'lat', in this order, and any number of layers to downscale- each as a new column
#' @param covar.ras= a raster or raster stack of high-resolution raster layers at resolution and extent of downscaling.  These layers will be used as co-variates and careful consideration should be given properly selecting these.  
#' @param ncores number of CPUs to use for interpolating - each core will be assigned a layer to interpolate. Best to first try with a few cores- this process can require a lot of memory
#' @param tps if tps=TRUE then the residuals of the best model will be interpolated using a thin-plate-spline and the final downscaled layer will be adjusted by this layer (this is what ANUSPLIN does)
#' @param smooth.outputs.only if smooth.outputs.only=TRUE, this removes Boosted Regressive Trees and Random Forests for the options of algorithms.  Both occasionally produce 'blocky' outputs.  I always recommend first using smooth.outputs.only=FALSE for the intial analyses. Then if you are unhappy with the visual appearence of the created layers, consider smooth.outputs.only=TRUE. Be aware that now your model performance might dramatically decline because those two algorithms have been excluded
#' @return  This script interpolates noisy multi-variate data through machine learning ensembling using six algorithms: boosted regression trees (BRT), neural networks (NN); generalized additive model (GAM), multivariate adaptive regression splines (MARS), support vector machines (SVM) and random forests (RF). This function evaluates (via k-fold cross validation, where k=10) a method’s ability to predict the input data and ensembles of all combinations of the six algorithms weighting each from 0-1 and evaluting fit.  The best model will have the lowest AICc (with the number of parameters in AICc calculation corresponding the number of models in ensemble).  After the best model is determined, the function will run the ensemble on the full dataset.  Then residuals will be calculated and interpolated using thin-plate-smoothing splines, which will secondarily correct the final ensemble model. This package is a free open-source machine learning analog to the expensive ANUSPLIN software. To output final R2 values, model weights, algorithm(s) used, and rasters for use in GIS; use the 'machisplin.write.geotiff' function.  To output residuals use 'machisplin.write.residuals' and to output model loadings use 'machispline.write.loadings'.
#' @export
#' @import terra
#' @import snow
#' @import snowfall
#' @importFrom fields Tps
#' @importFrom randomForest randomForest
#' @importFrom earth earth
#' @importFrom earth evimp
#' @importFrom kernlab ksvm
#' @importFrom mgcv gam
#' @importFrom NeuralNetTools garson
#' @importFrom nnet nnet
#' @importFrom optimx optimx
#' @importFrom breakDown broken
#' @export

#' @examples
#' ######## EXAMPLE 1 ########
#' library(MACHISPLIN)
#' library(terra)
#' 
#' # Import a csv as shapefile:
#' Mydata<-read.delim("sampling.csv", sep=",", h=T)
#' 
#' #load rasters to use as high resolution co-variates for downscaling
#' ALT = terra::rast("SRTM30m.tif")
#' SLOPE = terra::rast("ln_slope.tif")
#' ASPECT = terra::rast("aspect.tif")
#' GEOMORPH = terra::rast("geomorphons.tif")
#' TWI = terra::rast("TWI.tif")
#' 
#' # function input: raster brick of covarites
#' raster_stack<-terra::c(ALT,SLOPE,TWI,GEOMORPH, ASPECT)
#' 
#' #run an ensemble machine learning thin plate spline 
#' interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2, tps=FALSE)
#' 
#' #to plot results (change number to select different output raster)
#' plot(interp.rast[[1]]$final)
#' 
#' #to view residuals (change number to select different output raster)
#' interp.rast[[1]]$residuals
#' 
#' #to view model loadings
#' interp.rast[[1]]$var.imp
#' 
#' #to view other model parameters and other parameters
#' interp.rast[[1]]$summary
#' 
#' ######## EXAMPLE 2 ########
#' library(MACHISPLIN)
#' library(raster)
#' 
#' ##load spatial data with (coordinates named exactly as 'long' and 'lat') and any number of layers to downscale
#' data(sampling)
#' Mydata<-sampling
#' 
#' #load rasters to use as high resolution co-variates for downscaling
#' ALT = raster(system.file("extdata", "alt.tif", package="MACHISPLIN"))
#' SLOPE = raster(system.file("extdata", "slope.tif", package="MACHISPLIN"))
#' TWI = raster(system.file("extdata", "TWI.tif", package="MACHISPLIN"))
#' 
#' # function input: raster brick of covarites
#' raster_stack<-stack(ALT,SLOPE,TWI)
#' 
#' #run an ensemble machine learning thin plate spline 
#' interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2, smooth.outputs.only=TRUE)
#' 
#' machisplin.write.geotiff(mltps.in=interp.rast)
#' machisplin.write.residuals(mltps.in=interp.rast)
#' machisplin.write.loadings(mltps.in=interp.rast)

machisplin.mltps<-function(int.values, covar.ras, n.cores=1, tps=TRUE, smooth.outputs.only=FALSE, trouble=FALSE){
f <- list()
##TEMP OVERRIDE
n.cores=1
i.lyrs<-ncol(int.values)
# n of iterpolations - subtrack x and y
n.spln<-i.lyrs-2

# n of iterpolations-add two layers for lat and long
n.covars<-terra::nlyr(covar.ras)+2

## create a raster for lat and long and any other covariate, add to raster stack
## create rasters of LAT and LONG for downscalling, using ALT as template for dimensions
r<-covar.ras[[1]]
latitude_values <- terra::yFromCell(r, 1:terra::ncell(r))
LAT <- terra::setValues(r, latitude_values)
names(LAT) <- "LAT"
longitude_values <- terra::xFromCell(r, 1:terra::ncell(r))
LONG <- terra::setValues(r, longitude_values)
names(LONG) <- "LONG"

#########################################
#load raster for statistical downscaling (a raster brick saved in the grid format and countaining all topographic data)
rast.names<-names(covar.ras)
rast_stack<-c(covar.ras,LONG, LAT)
rast.names.all<-c(rast.names,"LONG","LAT")
n.names<-(1:terra::nlyr(rast_stack))
for(i in n.names){
names(rast_stack)[[i]]<- rast.names.all[i]}

#extract values to points from rasters
RAST_VAL<-data.frame(terra::extract(rast_stack, int.values[1:2]))
#str(RAST_VAL)

#merge sampled data to input
Mydata<-cbind(int.values,RAST_VAL)

#store full dataset rows
Mydata.FULL<-nrow(Mydata)
#remove NAs
Mydata<-Mydata[complete.cases(Mydata),]

#error term- report percent of rows with NA
if(nrow(Mydata)/Mydata.FULL<0.75){print(paste("Warning!",(Mydata.FULL-nrow(Mydata)), " points fell outside of input co-variate rasters (of",Mydata.FULL,"total input).  Consider using co-variates that match the full extent of the input data"))}

#Convert to spatial data frame ~ shapefile in R: SpatialPointsDataFrame (data"XY",data)
MyShp<-sp::SpatialPointsDataFrame(Mydata[,1:2],Mydata)

#specify coordinate system
terra::crs(MyShp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
 
#######################
#make list to loop over 
#the column 3 to 4 in temp_dat are the temperature variable 
#and the column 5 to xx are the topographic variables 
dat_tps<-list()
n.bios<-c(1:n.spln)# number of climate variables, here two
for (i in n.bios){
    dat_tps[[i]]<-MyShp[,c((i+2), c((i.lyrs+2):(n.covars+i.lyrs+1)))] #2 is b/c 3 and 4 col are temp data thus i+2=3 & 4 columns
  }
#str(dat_tps)

#rename each dataset within the list
name.dat<-c("resp",rast.names.all)
dat_tps<- lapply(seq(dat_tps), function(i) {
        y <- data.frame(dat_tps[[i]])
		yy=n.covars+1
		y <- y[1:yy]
        names(y) <- c(name.dat)
        return(y)
    })
#str(dat_tps)
dat_tps<<-dat_tps
#extract names for hi-res rasters
nm.spln.out<-c(3:(n.spln+2))
#out.names<-names(Mydata[,nm.spln.out])
out.names<-names(Mydata[nm.spln.out])


#extract names for model formula
xnam <- name.dat[2:(n.covars+1)]
mod.form<- as.formula(paste("resp ~ ", paste(xnam, collapse= "+")))

########################################################################
#########################################################################
#function to pass to snowfall for parallel computing 
if(n.cores>1){
	if(n.spln>1){i<-seq(1,n.spln)} else {i<-1}# length = number of climate variables
		myLapply_elevOut<-function(i){
			##################################################################################################
			############################# part 1 evaluate best ensemble of models ##############################
			##################################################################################################
            #perform K-fold cross val
			#l <- list()
			nfolds <- 10
			kfolds <- machisplin.kfold(dat_tps[[i]], nfolds)
						
			mfit.brt.full<-mfit.rf.full<-mfit.nn.full<-mfit.mars.full<-mfit.svm.full<-mfit.gam.full<-NULL
            
			for (v in 1:nfolds) {
			#train with 90% of data is <4000 or 10% if >4000 total rows
			if (nrow(dat_tps[[i]])>4000){
				train <- dat_tps[[i]][kfolds==v,]
				test <- dat_tps[[i]][kfolds!=v,]
			} else {train <- dat_tps[[i]][kfolds!=v,]
				test <- dat_tps[[i]][kfolds==v,]}
			
			##### create dataset for nueral networks (resp has to have values between 0-1)
			min.resp<-min(train$resp)
			max.resp<-max(train$resp)
			trainNN<-train
			if(min.resp<0){trainNN$resp<-trainNN$resp-min.resp}
			if(min.resp>=0){trainNN$resp<-trainNN$resp-min.resp}
			max2.resp<-max(trainNN$resp)
			trainNN$resp<-trainNN$resp/max2.resp
			
			####train models
            gc()
			mod.brt.tps.elev<- machispln.gbm.step(data=train, gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 25, learning.rate = 0.01, bag.fraction = 0.5, plot.main = FALSE)
		    mod.rf.tps.elev<- ranndomForest::randomForest(mod.form, data = train)
            mod.nn.tps.elev<-nnet::nnet(mod.form, data = trainNN, size=10, linout=TRUE, maxit=10000)
			mod.mars.tps.elev<-earth::earth(mod.form, data = train, nfold=10)
            mod.svm.tps.elev<- kernlab::ksvm(mod.form, data = train)
			mod.gam.tps.elev<- mgcv::gam(mod.form, data = train)
			
			
			#extract residuals and calculate residual sum of squares
			#BRT
			gc()
			pred.brt.obs <- predict(mod.brt.tps.elev, test, n.trees=mod.brt.tps.elev$gbm.call$best.trees, type="response")
            res.brt.elev<-test[,1]-pred.brt.obs
			if(is.null(mfit.brt.full)==FALSE){
			    mfit.brt.r2<-res.brt.elev
			    mfit.brt.full<-c(mfit.brt.full,mfit.brt.r2)}	
			if(is.null(mfit.brt.full)==TRUE){mfit.brt.full<-res.brt.elev}
	        	
			#RF
			gc()
            pred.rf.obs<-predict(mod.rf.tps.elev, test, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
			res.rf.elev<-test[,1]-pred.rf.obs
			if(is.null(mfit.rf.full)==FALSE){
			    mfit.rf.r2<-res.rf.elev
			    mfit.rf.full<-c(mfit.rf.full,mfit.rf.r2)}
			if(is.null(mfit.rf.full)==TRUE){mfit.rf.full<-res.rf.elev}
				
			#NN
			gc()
			pred.nn1<-predict(mod.nn.tps.elev, test)*max2.resp
			pred.nn2<-pred.nn1+min.resp
			res.nn.elev<-test[,1]-pred.nn2
			if(is.null(mfit.nn.full)==FALSE){
			    mfit.nn.r2<-res.nn.elev
			    mfit.nn.full<-c(mfit.nn.full,mfit.nn.r2)}
			if(is.null(mfit.nn.full)==TRUE){mfit.nn.full<-res.nn.elev}
            
			#MAR
			gc()
			pred.mars.test<-predict(mod.mars.tps.elev,test)
			res.mars.elev<-as.vector(test[,1]-pred.mars.test) 
			if(is.null(mfit.mars.full)==FALSE){
			    mfit.mars.r2<-res.mars.elev
			    mfit.mars.full<-c(mfit.mars.full,mfit.mars.r2)}
			if(is.null(mfit.mars.full)==TRUE){mfit.mars.full<-res.mars.elev}
            
			#SVM
			gc()
			pred.svm.test<-predict(mod.svm.tps.elev,test)
			res.svm.elev<-test[,1]-pred.svm.test 
			if(is.null(mfit.svm.full)==FALSE){
			    mfit.svm.r2<-res.svm.elev
			    mfit.svm.full<-c(mfit.svm.full,mfit.svm.r2)}
			if(is.null(mfit.svm.full)==TRUE){mfit.svm.full<-res.svm.elev}
            
			#GAM
			gc()
			pred.gam.test<-predict(mod.gam.tps.elev,test)
			res.gam.elev<-as.vector(test[,1]-pred.gam.test)
			if(is.null(mfit.gam.full)==FALSE){
			    mfit.gam.r2<-res.gam.elev
			    mfit.gam.full<-c(mfit.gam.full,mfit.gam.r2)}
			if(is.null(mfit.gam.full)==TRUE){mfit.gam.full<-res.gam.elev}
			}
			
			#evalute weighted ensemble
			#pick best weighted model
			if(smooth.outputs.only==FALSE){
				par<-c(0.5,0.5,0.5,0.5,0.5,0.5)
                
				machisplin.optimx.internal<- function(k1,k2,k3,k4,k5,k6,mfit.brt.full,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.rf.full,mfit.svm.full){
				fit<-sum((((mfit.brt.full*k1)/(k1+k2+k3+k4+k5+k6))+((mfit.gam.full*k2)/(k1+k2+k3+k4+k5+k6))+((mfit.nn.full*k3)/(k1+k2+k3+k4+k5+k6))+((mfit.mars.full*k4)/(k1+k2+k3+k4+k5+k6))+((mfit.rf.full*k5)/(k1+k2+k3+k4+k5+k6))+((mfit.svm.full*k6)/(k1+k2+k3+k4+k5+k6)))^2)
				return(fit)}
				
				OptX<-optimx::optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4], x[5], x[6],mfit.brt.full,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.rf.full,mfit.svm.full), lower=0, upper=1, method = "L-BFGS-B")
				
				OptX.mfit.wt<-OptX.mfit<-OptX.mfit.txt<-NULL
				OptX.mfit.wt.tot<-OptX$p1+OptX$p2+OptX$p3+OptX$p4+OptX$p5+OptX$p6
				OptX.cut<-0.05*OptX.mfit.wt.tot
				
					if(round(OptX$p1,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"b")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p1,2))}
						
				if(round(OptX$p2,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"g")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p2,2))}
	
				if(round(OptX$p3,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"n")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p3,2))}
						
				if(round(OptX$p4,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"m")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p4,2))}

				if(round(OptX$p5,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"r")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p5,2))}
					
				if(round(OptX$p6,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"v")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p6,2))}

			
			}
			if(smooth.outputs.only==TRUE){
				par<-c(0.5,0.5,0.5,0.5)
				
				machisplin.optimx.internal<- function(k2,k3,k4,k6,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.svm.full){
				fit<-sum((((mfit.gam.full*k2)/(k2+k3+k4+k6))+((mfit.nn.full*k3)/(k2+k3+k4+k6))+((mfit.mars.full*k4)/(k2+k3+k4+k6))+((mfit.svm.full*k6)/(k2+k3+k4+k6)))^2)
				return(fit)}
				
				OptX<-optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4],mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.svm.full), lower=0, upper=1, method = "L-BFGS-B")
				
				OptX.mfit.wt<-OptX.mfit<-OptX.mfit.txt<-NULL
				OptX.mfit.wt.tot<-OptX$p1+OptX$p2+OptX$p3+OptX$p4
				OptX.cut<-0.05*OptX.mfit.wt.tot
				if(round(OptX$p1,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"g")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p1,2))}
					
				if(round(OptX$p2,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"n")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p2,2))}
				
				if(round(OptX$p3,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"m")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p3,2))}
					
				if(round(OptX$p4,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"v")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p4,2))}
			}
			#l$ensemble.models<-OptX.mfit
			#l$ensemble.weights<-OptX.mfit.wt
			f.mod<-OptX.mfit

			#get number of models in best
            ku<-nchar(f.mod)
			
			#create vector of models
			k=c(1:ku)
			iter.mod<-0
		    #split out letters of each model
			mods.run<-unlist(strsplit(f.mod,""))
			OptX.per.tot<-sum(OptX.mfit.wt)
			
			for(k in mods.run){	
       		    iter.mod<-iter.mod+1
				if(k=="b"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="g"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="n"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="m"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="r"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="v"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
			     }
				if(OptX.mfit.txt==1){OptX.mfit.txt="none"}
			##################################################################################################
			######################################## part 2 run best model(s)  #################################
			##################################################################################################
            
			pred.elev<- NULL
			res.FINAL<-NULL
			iter.mod<-0
			for(k in mods.run){	
			  iter.mod<-iter.mod+1
   			  if (k=="n"){
			    ##### create dataset for nueral networks (resp has to have values between 0-1)
			    gc()
				trainNN.f<-dat_tps[[i]]
			    min.resp.f<-min(trainNN.f$resp)
			    max.resp.f<-max(trainNN.f$resp)
			    trainNN.f$resp<-trainNN.f$resp-min.resp.f
			    max2.resp.f<-max(trainNN.f$resp)
			    trainNN.f$resp<-trainNN.f$resp/max2.resp.f
                
			    mod.run="NN"
			    #run model with all points
			    mod.nn.tps.FINAL<-nnet::nnet(mod.form, data = trainNN.f, size=10, linout=TRUE, maxit=10000)
			    #store variable importance
                l$var.imp$nn<-NeuralNetTools::garson(mod.nn.tps.FINAL, bar_plot=F)
				
			    #create raster pred
			    if(is.null(pred.elev)==FALSE){pred.elev.nn<-predict(rast_stack, mod.nn.tps.FINAL)
			    pred.nn<-pred.elev.nn*max2.resp.f
		        pred.elev.2<-pred.nn+min.resp.f
			    pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
			    if(is.null(pred.elev)==TRUE){pred.elev.nn<-predict(rast_stack, mod.nn.tps.FINAL)
			    pred.nn<-pred.elev.nn*max2.resp.f
			    pred.elev.2<-pred.nn+min.resp.f
                pred.elev<-(pred.elev.2*OptX.mfit.wt[[iter.mod]])}
			    #get residuals
			    pred.resid.nn<-predict(mod.nn.tps.FINAL,trainNN.f)
				pred.resid.nn<-pred.resid.nn*max2.resp.f
				pred.resid.nn<-pred.resid.nn+min.resp.f
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-dat_tps[[i]]$resp-pred.resid.nn
				res.FINAL<-res.FINAL+(res.FINAL.2*OptX.mfit.wt[[iter.mod]])}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(dat_tps[[i]]$resp-pred.resid.nn)*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
             			
  	          ##final models = boosted regression tree	
			  if (k=="b"){
				gc()
				mod.run="BRT"
				#run model with all points
				mod.brt.tps.FINAL<- machisplin.gbm.step(data=dat_tps[[i]], gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5, plot.main = FALSE)
				#store variable importance
				l$var.imp$brt<-mod.brt.tps.FINAL$contributions
				#create raster prediction
				if(is.null(pred.elev)==FALSE){pred.elev.2<- predict(rast_stack, mod.brt.tps.FINAL, n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")
				pred.elev<-pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]])}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.brt.tps.FINAL, n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")*OptX.mfit.wt[[iter.mod]]}
				#predict at all train sites
				pred.brt.obs <- predict(mod.brt.tps.FINAL, dat_tps[[i]], n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")
				#get residuals
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-(dat_tps[[i]]$resp-pred.brt.obs)*OptX.mfit.wt[[iter.mod]]
				res.FINAL<-res.FINAL+res.FINAL.2}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(dat_tps[[i]]$resp-pred.brt.obs)*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
			  
  		      ##final models = randomForest	
         	  #if (mfit.rf<mfit.brt & mfit.rf<mfit.lm & mfit.rf=<mfit.mars){
			  if (k=="r"){
				gc()
				mod.run="RF"
				#run model with all points
				mod.rf.tps.FINAL<- randomForest::randomForest(mod.form, data = dat_tps[[i]],importance = TRUE)
				#store variable importance
				l$var.imp$rf<-mod.rf.tps.FINAL$importance
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.rf.tps.FINAL, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
				pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.rf.tps.FINAL, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				pred.rf.obs<-predict(mod.rf.tps.FINAL, dat_tps[[i]], type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-(dat_tps[[i]]$resp-pred.rf.obs)*OptX.mfit.wt[[iter.mod]]
				res.FINAL<-(res.FINAL+res.FINAL.2)}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(dat_tps[[i]]$resp-pred.rf.obs)*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
                
			  ##final models = MARS	
			  if (k=="m"){
				gc()
				mod.run="MARS"
				#run model with all points
				mod.MARS.tps.FINAL<-earth::earth(mod.form,  data = dat_tps[[i]], nfold=10)
				#store variable importance
				l$var.imp$mars<-earth::evimp(mod.MARS.tps.FINAL)
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.MARS.tps.FINAL)
				pred.elev<-pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]])}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.MARS.tps.FINAL)*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-resid(mod.MARS.tps.FINAL, warn=FALSE)
				res.FINAL<-(res.FINAL+(as.vector(res.FINAL.2)*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(as.vector(resid(mod.MARS.tps.FINAL, warn=FALSE)))*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
              
   		      ##final models = SVM
			  if (k=="v"){
				gc()
				mod.run="SVM"
				#run model with all points
				mod.SVM.tps.FINAL<-kernlab::ksvm(mod.form, data=dat_tps[[i]])
				#store variable importance
				nvals.in<-nrow(dat_tps[[i]])
				sampleSVM<-dat_tps[[i]]
				if(nrow(sampleSVM)>200){sampleSVM<- sampleSVM[sample(nrow(sampleSVM),200),]}
				nvals.in<-nrow(sampleSVM)
				SVM.imp<-NULL
				for(sub.svm in 1:nvals.in){
				nob <- sampleSVM[sub.svm,]
				nob <- nob[1:(1+n.covars)]
				set.seed(1313)	
				explain_z<- broken(mod.SVM.tps.FINAL,  new_observation = nob, data = sampleSVM, predict.function = predict, baseline = "intercept", direction = "up") 
				if(is.null(SVM.imp)==TRUE){SVM.imp<-abs(explain_z$contribution)}
				if(is.null(SVM.imp)==FALSE){SVM.imp<-abs(SVM.imp)+abs(explain_z$contribution)}
				}
				SVM.f<-(SVM.imp/nvals.in)
				SVM.f<-as.matrix(SVM.f,ncol=1)
				name.SVM.imp<-gsub(" =.*","",explain_z$variable)
				rownames(SVM.f)<-name.SVM.imp
				colnames(SVM.f)<-("contributions to SVM")
				l$var.imp$SVM<-SVM.f
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.SVM.tps.FINAL)
				pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.SVM.tps.FINAL)*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				pred.SVM.obs<-predict(mod.SVM.tps.FINAL, dat_tps[[i]])
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-as.vector(dat_tps[[i]]$resp-pred.SVM.obs)
				res.FINAL<-(res.FINAL+(res.FINAL.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(as.vector(dat_tps[[i]]$resp-pred.SVM.obs))*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
              
   		      ##final models = GAM
			  if (k=="g"){
				gc()
				mod.run="GAM"
				#run model with all points
				mod.GAM.tps.FINAL<-mgcv::gam(mod.form, data=dat_tps[[i]])
				#store variable importance
				l$var.imp$gam<-mod.GAM.tps.FINAL$coefficients
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.GAM.tps.FINAL)
				pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(pred.elev)==TRUE){pred.elev<-(predict(rast_stack, mod.GAM.tps.FINAL))*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				pred.GAM.obs<-predict(mod.GAM.tps.FINAL, dat_tps[[i]])
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-as.vector(dat_tps[[i]]$resp-pred.GAM.obs)
				res.FINAL<-(res.FINAL+(res.FINAL.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(as.vector(dat_tps[[i]]$resp-pred.GAM.obs))*OptX.mfit.wt[[iter.mod]]}
	   		    gc()
				}
			}
			#wrap up analysis and get final layer
			#divide models by number ensembled
            pred.elev<-(pred.elev/OptX.mfit.wt.tot)
			res.FINAL<-(res.FINAL/OptX.mfit.wt.tot)
			
			#calculate Sum of squares and resisdual sum of squares
			rss.m <- sum((res.FINAL) ^ 2)
			tss <- sum((dat_tps[[i]]$resp - mean(dat_tps[[i]]$resp)) ^ 2)
			l$residuals<- cbind((res.FINAL),dat_tps[[i]]$LONG,dat_tps[[i]]$LAT)
			colnames(l$residuals)<-c("residuals","long","lat")

			rsq.model <- 1 - (rss.m/tss) #r squared
			gc()
			##################################################################################################
			###################################### part 3 spline residuals ###################################
			##################################################################################################
			#DETERMINE IF TILING IS NEEDED
			if(tps==TRUE) {
				#specify total extent
				totalExt<-extent(rast_stack)
				#specify n rows 
				nr.in<-nrow(rast_stack)
				#specify n cols 
				nc.in<-ncol(rast_stack)
				#specify n row blocks
				nrB<-nr.in/3000
				nRx<-ceiling(nrB)
				#specify n col blocks
				ncB<-nc.in/3000
				nCx<-ceiling(ncB)
				#specify lat distance
				longDist<-((totalExt[2]-totalExt[1])/nCx)
				#specify long distance
				latDist<-((totalExt[4]-totalExt[3])/nRx)
				
				#specify sampling grid for TPS
				m2<-0;m1<-0;new.df<-NULL;new.df2<-NULL
				
				for (j in 1:nRx){
				for (h in 1:nCx){
						m1<-m1+1
						new.df[[m1]]<-c((totalExt[1]+((longDist*(h-1))-(longDist*0.2))),(totalExt[1]+((longDist*h)+(longDist*0.2))),(totalExt[3]+((latDist*(j-1)))-(latDist*0.2)),(totalExt[3]+((latDist*j))+(latDist*0.2)))
					}
				}
				#specify sampling grid for mosaic
				for (j in 1:nRx){
				   for (h in 1:nCx){
						m2<-m2+1
						new.df2[[m2]]<-c((totalExt[1]+((longDist*(h-1))-(longDist*0.025))),(totalExt[1]+((longDist*h)+(longDist*0.025))),(totalExt[3]+((latDist*(j-1)))-(latDist*0.025)),(totalExt[3]+((latDist*j))+(latDist*0.025)))
				}}
				
				#specify full cordinates for spatial subsampling
				#perform TPS on each grid and then mosaic
				if(nRx*nCx>1){
					Full.cords<-dat_tps[[1]][,c(n.covars,n.covars+1)]
					pred_TPS_elev<-NULL
					for (h in 1:(nRx*nCx)){
						gc()
						#clip raster brick
						b <- as(extent(new.df[[h]][1], new.df[[h]][2], new.df[[h]][3], new.df[[h]][4]), 'SpatialPolygons')
						d <- as(extent(new.df2[[h]][1], new.df2[[h]][2], new.df2[[h]][3], new.df2[[h]][4]), 'SpatialPolygons')
						crs(b) <- crs(rast_stack)
						crs(d) <- crs(rast_stack)
						rb <- crop(rast_stack, b)
						#sub sample residuals
						RAST_TPS<-data.frame(extract(rb[[1]], Full.cords))
						#str(RAST_VAL)
						#merge sampled data to input
						MyTPSdata<-cbind(res.FINAL,Full.cords,RAST_TPS)
						#remove NAs
						MyTPSdata<-MyTPSdata[complete.cases(MyTPSdata),][1:3]
						print(paste("Number of data points in grid:",nrow(MyTPSdata)))
						#fit thin plate spline of residuals
						if(nrow(MyTPSdata)<10){
							print("Sample grid mostly empty: grid area likely is water or a region with no input data. This grid will not be part of TPS. This is expected in many situations and part of the error correction process.")
							rbT<-rb[[1]]
							r<-raster(rbT)
							r[is.na(r[])] <- 0
							pred_TPS_elev<-r
							TPS_name<-paste0("TILE_",h)
							names(pred_TPS_elev)<- TPS_name
							pred_TPS_elev <- crop(pred_TPS_elev, d)
							pred_TPS_elev <- extend(pred_TPS_elev,rast_stack)
							plot(pred_TPS_elev)
							Sys.sleep(0)}
						else {mod.tps.elev<-Tps(MyTPSdata[2:3], MyTPSdata[1])#columns of lat and long
							#use TPS to interpolate residuals
							TPS_name<-paste0("TILE_",h)
							#use TPS to interpolate residuals
							pred_TPS_elev<-interpolate(rb, mod.tps.elev)
							names(pred_TPS_elev)<- TPS_name
							pred_TPS_elev <- crop(pred_TPS_elev, d)
							pred_TPS_elev <- extend(pred_TPS_elev,rast_stack)
							plot(pred_TPS_elev)
							Sys.sleep(0)}
						if (h>1){raster_sTPS<-stack(raster_sTPS,pred_TPS_elev)}
						else {raster_sTPS<-pred_TPS_elev}
						gc()
						}}
				if(nRx*nCx>1){
					rast.list<-as.list(raster_sTPS)
					rast.list$fun <- mean
					rast.mosaic <- do.call(mosaic,rast.list)}
				if(nRx*nCx==1){
					#fit thin plate spline of residuals
					mod.tps.elev<-Tps(dat_tps[[i]][,c(n.covars,n.covars+1)], res.FINAL)#columns of lat and long
					#use TPS to interpolate residuals
					rast.mosaic<-interpolate(rast_stack, mod.tps.elev)}
				##################################################################################################
				############################## part 4 feather seams of spline tiles ##############################
				##################################################################################################
				#feather right
				feath.ras.TPS<-NULL
				
				for (j in 1:nRx){
				   for (h in 1:nCx){
						gc()
						v<-(h+((j*nCx)-nCx))
						if (h<nCx){
						AAA<-raster_sTPS[[v]]+raster_sTPS[[v+1]]
						BBB<-as.data.frame(rasterToPoints(AAA))
						#Convert to spatial data frame ~ shapefile in R: SpatialPointsDataFrame (data"XY",data)
						MyShp<-SpatialPointsDataFrame(BBB[,1:2],BBB)
						
						#specify coordinate system
						crs(MyShp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
						#
						blend.ext<-extent(MyShp)
						#crop
						bR1 <- crop(raster_sTPS[[v]], blend.ext)
						bR2 <- crop(raster_sTPS[[v+1]], blend.ext)
						
						LONG.B1 <- bR1
						xy <- coordinates(bR1)
						LONG.B1[] <- xy[, 1]
						LONG.B2 <- bR2
						xy <- coordinates(bR2)
						LONG.B2[] <- xy[, 1]
						
						delta.L1<-(LONG.B1@data@max)-(LONG.B1@data@min)
						stD1<-(LONG.B1-LONG.B1@data@min)/delta.L1
						stD1<-1-stD1
						
						delta.L2<-(LONG.B2@data@max)-(LONG.B2@data@min)
						stD2<-(LONG.B2-LONG.B2@data@min)/delta.L2
						
						we.R.d1<-stD1*bR1
						we.R.d2<-stD2*bR2
						feath.ras<-we.R.d2+we.R.d1
						feath.ras<- extend(feath.ras,rast_stack)
						if(is.null(feath.ras.TPS)==TRUE){feath.ras.TPS<-feath.ras}else{feath.ras.TPS<-stack(feath.ras,feath.ras.TPS)}
						gc()
						}
						#####################################################
						#####################################################
						#feather up
						if (j<nRx){
						gc()
						
						AAA<-raster_sTPS[[v]]+raster_sTPS[[v+nCx]]
						BBB<-as.data.frame(rasterToPoints(AAA))
						#Convert to spatial data frame ~ shapefile in R: SpatialPointsDataFrame (data"XY",data)
						MyShp<-SpatialPointsDataFrame(BBB[,1:2],BBB)
						
						#specify coordinate system
						crs(MyShp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
						
						#
						blend.ext<-extent(MyShp)
						#crop
						bR1 <- crop(raster_sTPS[[v]], blend.ext)
						bR2 <- crop(raster_sTPS[[v+nCx]], blend.ext)
						LAT.B1 <- bR1
						xy <- coordinates(bR1)
						LAT.B1[] <- xy[, 2]
						
						LAT.B2<- bR2
						xy <- coordinates(bR2)
						LAT.B2[] <- xy[, 2]
						
						delta.L1<-(LAT.B1@data@max)-(LAT.B1@data@min)
						stD1<-(LAT.B1-LAT.B1@data@min)/delta.L1
						stD1<-1-stD1
						
						delta.L2<-(LAT.B2@data@max)-(LAT.B2@data@min)
						stD2<-(LAT.B2-LAT.B2@data@min)/delta.L2
						
						we.R.d1<-stD1*bR1
						we.R.d2<-stD2*bR2
						feath.ras<-we.R.d2+we.R.d1
						feath.ras<- extend(feath.ras,rast_stack)
						if(is.null(feath.ras.TPS)==TRUE){feath.ras.TPS<-feath.ras}else{feath.ras.TPS<-stack(feath.ras,feath.ras.TPS)}
						gc()
						}}}
				if(nRx*nCx>2){
					rast.list<-as.list(feath.ras.TPS)
						rast.list$fun <- mean
								rast.mosaic.TPS <- do.call(mosaic,rast.list)
					#r.TPS<-rast.mosaic.TPS
					#r.rastmo.min<-rast.mosaic@data@min
					#r.TPS[!is.na(r.TPS[])] <- r.rastmo.min
					final.TPS<-merge(rast.mosaic.TPS,rast.mosaic)}
				if(nRx*nCx==2){
					final.TPS<-merge(feath.ras.TPS,rast.mosaic)}
				if(nRx*nCx==1){
					final.TPS<-rast.mosaic}
			}
			##################################################################################################
			################################# part 5 finalize results ##################################
			##################################################################################################
			
			if(tps==TRUE) {
				gc()
				#calculate pred at normal scale, suming up kriging pred+res
				pred.elev.i<- brick(pred.elev, final.TPS)
				pred.elev.i.calc<- calc(pred.elev.i, fun=sum)
				
				#extract final values to input points from final raster
				f.actual<-extract(pred.elev.i.calc, dat_tps[[i]][,n.covars:(n.covars+1)])#lat and long input
            
				#calculate Sum of squares and resisdual sum of squares
				rss.final <- sum((dat_tps[[i]]$resp - f.actual) ^ 2)
				l$resid_TPS<- cbind((dat_tps[[i]]$resp - f.actual),dat_tps[[i]]$LONG,dat_tps[[i]]$LAT)
				colnames(l$resid_TPS)<-c("residuals","long","lat")
				rsq.final <- 1 - (rss.final/tss) #r squared
				
				#out values
				sum.val<-data.frame(out.names[i],f.mod,OptX.mfit.txt,rsq.model,rsq.final)
				colnames(sum.val)<-c("layer","best model(s):","ensemble weights:","r2 ensemble:","r2 final:")
				if(i==1){l$n.layers<-length(out.names)}
				l$summary<-sum.val
				#save TPS if R2 is better, else save only model
				if(rsq.final>rsq.model){
					names(pred.elev.i.calc)<- out.names[i]
					l$final<-pred.elev.i.calc
					} else {
					names(pred.elev)<- out.names[i]
					l$final<-pred.elev
					}
			}
			if(tps==FALSE) {
				gc()
				#out values
				sum.val<-data.frame(out.names[i],f.mod,OptX.mfit.txt,rsq.model)
				colnames(sum.val)<-c("layer","best model(s):","ensemble weights:","r2 ensemble:")
				l$n.layers<-length(out.names)
				l$summary<-sum.val
				#save TPS if R2 is better, else save only model
				names(pred.elev)<- out.names[i]
				l$final<-pred.elev
				}
		l$var.imp
		return(l)
    }
 
          
	#Initiate snowfall
	sfInit(parallel=TRUE, cpus=n.cores)
	## Export packages
	sfLibrary('raster', warn.conflicts=FALSE, character.only=TRUE)
	#sfLibrary('gstat', warn.conflicts=FALSE, character.only=TRUE)
	#sfLibrary('nlme', warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('fields', warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('dismo', warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('randomForest', warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('earth',  warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('kernlab',  warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('mgcv',  warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('optimx',  warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('nnet',  warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('NeuralNetTools',  warn.conflicts=FALSE, character.only=TRUE)
	sfLibrary('breakDown',  warn.conflicts=FALSE, character.only=TRUE)
	
	## Export variables
	sfExport('dat_tps')
	sfExport('rast_stack')
	sfExport('i')
	sfExport('mod.form')
	sfExport('n.covars')
	sfExport('out.names')
	sfExport('tps')
	## Do the run
	mySFelevOut <- sfLapply(i, myLapply_elevOut)
	## stop snowfall
	#sfStop(nostop=FALSE)
	sfStop()
	f<-mySFelevOut
	return(f)
}
if(n.cores==1){
    sink('MachiSplin.LOG.txt', append=FALSE, split=TRUE)
	omega <- as.list(rep(NA, n.spln))
	if(n.spln>1){iter.clim<-seq(1,n.spln)} else {iter.clim<-1}# length = number of climate variables
		for(i in iter.clim){
			##################################################################################################
			############################# part 1 evaluate best ensemble of models ##############################
			##################################################################################################
            #perform K-fold cross val
			print("###########################################################################################")
			print("################################   STARTING  DOWNSCALING   ################################")
			print("###########################################################################################")
			print("##################################### STEP 1 (of 3) #######################################")
			print("###########################################################################################")
			print("################### Running Intial Models to Determine Best Algorithms ####################")
			print("###########################################################################################")
			print(paste("###########################  ", (out.names[i]), "  ######################################"))
			print("###########################################################################################")
			print("###########################################################################################")
			print(paste("Setting up datasets for k-fold crossvalidation: ",(out.names[i])))
			print(paste("                                         ",Sys.time()))
			nfolds <- 10
			kfolds <- machisplin.kfold(dat_tps[[i]], nfolds)
			l <- list()			
			mfit.brt.full<-mfit.rf.full<-mfit.nn.full<-mfit.mars.full<-mfit.svm.full<-mfit.gam.full<-NULL
			
			for (v in 1:nfolds) {
			print(paste("*******Iteration",(v),"(of 10)*************************************************"))
			#train with 90% of data is <4000 or 10% if >4000 total rows
			if (nrow(dat_tps[[i]])>4000){
				train <- dat_tps[[i]][kfolds==v,]
				test <- dat_tps[[i]][kfolds!=v,]
			} else {train <- dat_tps[[i]][kfolds!=v,]
				test <- dat_tps[[i]][kfolds==v,]}
			
			##### create dataset for nueral networks (resp has to have values between 0-1)
			min.resp<-min(train$resp)
			max.resp<-max(train$resp)
			trainNN<-train
			if(min.resp<0){trainNN$resp<-trainNN$resp-min.resp}
			if(min.resp>=0){trainNN$resp<-trainNN$resp-min.resp}
			max2.resp<-max(trainNN$resp)
			trainNN$resp<-trainNN$resp/max2.resp
			
			####train models
			print(paste("Running k-fold cross-validation for ",(out.names[i]),": Boosted Regresion Trees: Iteration",(v),"(of 10)"))
            print(paste("                                         ",Sys.time()))
			gc()
			mod.brt.tps.elev<- machisplin.gbm.step(data=train, gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 25, learning.rate = 0.01, bag.fraction = 0.5, plot.main = FALSE, silent = TRUE)
		    mod.rf.tps.elev<- randomForest::randomForest(mod.form, data = train)
            mod.nn.tps.elev<-nnet::nnet(mod.form, data = trainNN, size=10, linout=TRUE, maxit=10000)
			mod.mars.tps.elev<-earth::earth(mod.form, data = train, nfold=10)
            mod.svm.tps.elev<- kernlab::ksvm(mod.form, data = train)
			mod.gam.tps.elev<- mgcv::gam(mod.form, data = train)
			
			
			#extract residuals and calculate residual sum of squares
			#BRT
			gc()
			pred.brt.obs <- terra::predict(mod.brt.tps.elev, test, n.trees=mod.brt.tps.elev$gbm.call$best.trees, type="response")
            res.brt.elev<-test[,1]-pred.brt.obs
			if(is.null(mfit.brt.full)==FALSE){
			    mfit.brt.r2<-res.brt.elev
			    mfit.brt.full<-c(mfit.brt.full,mfit.brt.r2)}	
			if(is.null(mfit.brt.full)==TRUE){mfit.brt.full<-res.brt.elev}
	        	
			#RF
			gc()
			print(paste("Running k-fold cross-validation for ",(out.names[i]),": Random Forests: Iteration",(v),"(of 10)"))
			print(paste("                                         ",Sys.time()))
            pred.rf.obs<-terra::predict(mod.rf.tps.elev, test, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
			res.rf.elev<-test[,1]-pred.rf.obs
			if(is.null(mfit.rf.full)==FALSE){
			    mfit.rf.r2<-res.rf.elev
			    mfit.rf.full<-c(mfit.rf.full,mfit.rf.r2)}
			if(is.null(mfit.rf.full)==TRUE){mfit.rf.full<-res.rf.elev}
				
			#NN
			print(paste("Running k-fold cross-validation for ",(out.names[i]),": Neural Networks: Iteration",(v),"(of 10)"))
			print(paste("                                         ",Sys.time()))
			gc()
			pred.nn1<-terra::predict(mod.nn.tps.elev, test)*max2.resp
			pred.nn2<-pred.nn1+min.resp
			res.nn.elev<-test[,1]-pred.nn2
			if(is.null(mfit.nn.full)==FALSE){
			    mfit.nn.r2<-res.nn.elev
			    mfit.nn.full<-c(mfit.nn.full,mfit.nn.r2)}
			if(is.null(mfit.nn.full)==TRUE){mfit.nn.full<-res.nn.elev}
            
			#MAR
			print(paste("Running k-fold cross-validation for ",(out.names[i]),": Multivariate Adaptive Regression Splines: Iteration",(v),"(of 10)"))
			print(paste("                                         ",Sys.time()))
			gc()
			pred.mars.test<-terra::predict(mod.mars.tps.elev,test)
			res.mars.elev<-as.vector(test[,1]-pred.mars.test) 
			if(is.null(mfit.mars.full)==FALSE){
			    mfit.mars.r2<-res.mars.elev
			    mfit.mars.full<-c(mfit.mars.full,mfit.mars.r2)}
			if(is.null(mfit.mars.full)==TRUE){mfit.mars.full<-res.mars.elev}
            
			#SVM
			gc()
			print(paste("Running k-fold cross-validation for ",(out.names[i]),": Support Vector Machines: Iteration",(v),"(of 10)"))
			print(paste("                                         ",Sys.time()))
			pred.svm.test<-terra::predict(mod.svm.tps.elev,test)
			res.svm.elev<-test[,1]-pred.svm.test 
			if(is.null(mfit.svm.full)==FALSE){
			    mfit.svm.r2<-res.svm.elev
			    mfit.svm.full<-c(mfit.svm.full,mfit.svm.r2)}
			if(is.null(mfit.svm.full)==TRUE){mfit.svm.full<-res.svm.elev}
            
			#GAM
			gc()
			print(paste("Running k-fold cross-validation for ",(out.names[i]),": Generalized Additive Models: Iteration",(v),"(of 10)"))
			print(paste("                                         ",Sys.time()))
			pred.gam.test<-terra::predict(mod.gam.tps.elev,test)
			res.gam.elev<-as.vector(test[,1]-pred.gam.test)
			if(is.null(mfit.gam.full)==FALSE){
			    mfit.gam.r2<-res.gam.elev
			    mfit.gam.full<-c(mfit.gam.full,mfit.gam.r2)}
			if(is.null(mfit.gam.full)==TRUE){mfit.gam.full<-res.gam.elev}
			}
			
			#evalute weighted ensemble
			#pick best weighted model
			print("Finished Cross-Validation of Six Algorithms")
			print("Starting Ensemble Modeling")
			if(smooth.outputs.only==FALSE){
				par<-c(0.5,0.5,0.5,0.5,0.5,0.5)
                
				machisplin.optimx.internal<- function(k1,k2,k3,k4,k5,k6,mfit.brt.full,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.rf.full,mfit.svm.full){
				fit<-sum((((mfit.brt.full*k1)/(k1+k2+k3+k4+k5+k6))+((mfit.gam.full*k2)/(k1+k2+k3+k4+k5+k6))+((mfit.nn.full*k3)/(k1+k2+k3+k4+k5+k6))+((mfit.mars.full*k4)/(k1+k2+k3+k4+k5+k6))+((mfit.rf.full*k5)/(k1+k2+k3+k4+k5+k6))+((mfit.svm.full*k6)/(k1+k2+k3+k4+k5+k6)))^2)
				return(fit)}
				
				OptX<-optimx::optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4], x[5], x[6],mfit.brt.full,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.rf.full,mfit.svm.full), lower=0, upper=1, method = "L-BFGS-B")
				print("Best Ensemble Model Determined: considered all algorithms")

				OptX.mfit.wt<-OptX.mfit<-OptX.mfit.txt<-NULL
				OptX.mfit.wt.tot<-OptX$p1+OptX$p2+OptX$p3+OptX$p4+OptX$p5+OptX$p6
				OptX.cut<-0.05*OptX.mfit.wt.tot
				
					if(round(OptX$p1,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"b")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p1,2))}
						
				if(round(OptX$p2,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"g")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p2,2))}
	
				if(round(OptX$p3,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"n")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p3,2))}
						
				if(round(OptX$p4,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"m")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p4,2))}

				if(round(OptX$p5,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"r")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p5,2))}
					
				if(round(OptX$p6,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"v")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p6,2))}

			
			}
			if(smooth.outputs.only==TRUE){
				par<-c(0.5,0.5,0.5,0.5)
				
				machisplin.optimx.internal<- function(k2,k3,k4,k6,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.svm.full){
				fit<-sum((((mfit.gam.full*k2)/(k2+k3+k4+k6))+((mfit.nn.full*k3)/(k2+k3+k4+k6))+((mfit.mars.full*k4)/(k2+k3+k4+k6))+((mfit.svm.full*k6)/(k2+k3+k4+k6)))^2)
				return(fit)}
				print("Best Ensemble Model Determined: considered smooth outputs only")
				OptX<-optimx::optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4],mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.svm.full), lower=0, upper=1, method = "L-BFGS-B")
				
				OptX.mfit.wt<-OptX.mfit<-OptX.mfit.txt<-NULL
				OptX.mfit.wt.tot<-OptX$p1+OptX$p2+OptX$p3+OptX$p4
				OptX.cut<-0.05*OptX.mfit.wt.tot
				if(round(OptX$p1,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"g")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p1,2))}
					
				if(round(OptX$p2,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"n")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p2,2))}
				
				if(round(OptX$p3,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"m")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p3,2))}
					
				if(round(OptX$p4,2)>OptX.cut){
						OptX.mfit<-paste0(OptX.mfit,"v")
						OptX.mfit.wt<-c(OptX.mfit.wt, round(OptX$p4,2))}
			}
			#l$ensemble.models<-OptX.mfit
			#l$ensemble.weights<-OptX.mfit.wt
			f.mod<-OptX.mfit

			#get number of models in best
            ku<-nchar(f.mod)
			
			#create vector of models
			k=c(1:ku)
			iter.mod<-0
		    #split out letters of each model
			mods.run<-unlist(strsplit(f.mod,""))
			OptX.per.tot<-sum(OptX.mfit.wt)
			print("Setting up parameters to run ensemble models")
			for(k in mods.run){	
       		    iter.mod<-iter.mod+1
				if(k=="b"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="g"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="n"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="m"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="r"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
				if(k=="v"){
					if(is.null(OptX.mfit.txt)==FALSE){OptX.mfit.txt<-paste0(OptX.mfit.txt,":",round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1))}
					if(is.null(OptX.mfit.txt)==TRUE){OptX.mfit.txt<-round(((OptX.mfit.wt[[iter.mod]]/OptX.per.tot)*100),1)}}
			     }
				if(OptX.mfit.txt==1){OptX.mfit.txt="none"}
			##################################################################################################
			######################################## part 2 run best model(s)  #################################
			print("###########################################################################################")
			print("###########################################################################################")
			print("##################################### STEP 2 (of 3) #######################################")
			print("###########################################################################################")
			print(paste("############## Running Final Models for Ensemble of ",(ku)," Algorithms ###############"))
			print("###########################################################################################")
			print(paste("###########################  ", (out.names[i]), "  ######################################"))
			print("###########################################################################################")
			print("###########################################################################################")
			print(paste("                                         ",Sys.time()))
			pred.elev<- NULL
			res.FINAL<-NULL
			iter.mod<-0
			##TEMP TROUBLESHOOTING PARAM
			if (trouble==TRUE){mods.run<-"b"}
			for(k in mods.run){	
			  iter.mod<-iter.mod+1
   			  if (k=="n"){
			  	print(paste("Final modeling of ",(out.names[i]),": Nueral Networks"))
				print(paste("                                         ",Sys.time()))
			    ##### create dataset for nueral networks (resp has to have values between 0-1)
			    gc()
				trainNN.f<-dat_tps[[i]]
			    min.resp.f<-min(trainNN.f$resp)
			    max.resp.f<-max(trainNN.f$resp)
			    trainNN.f$resp<-trainNN.f$resp-min.resp.f
			    max2.resp.f<-max(trainNN.f$resp)
			    trainNN.f$resp<-trainNN.f$resp/max2.resp.f
                
			    mod.run="NN"
			    #run model with all points
			    mod.nn.tps.FINAL<-nnet::nnet(mod.form, data = trainNN.f, size=10, linout=TRUE, maxit=10000)
			    #store variable importance
                l$var.imp$nn<-NeuralNetTools::garson(mod.nn.tps.FINAL, bar_plot=F)
				
			    #create raster pred
			    if(is.null(pred.elev)==FALSE){pred.elev.nn<-terra::predict(rast_stack, mod.nn.tps.FINAL)
			    pred.nn<-pred.elev.nn*max2.resp.f
		        pred.elev.2<-pred.nn+min.resp.f
			    pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
			    if(is.null(pred.elev)==TRUE){pred.elev.nn<-terra::predict(rast_stack, mod.nn.tps.FINAL)
			    pred.nn<-pred.elev.nn*max2.resp.f
			    pred.elev.2<-pred.nn+min.resp.f
                pred.elev<-(pred.elev.2*OptX.mfit.wt[[iter.mod]])}
			    #get residuals
			    pred.resid.nn<-terra::predict(mod.nn.tps.FINAL,trainNN.f)
				pred.resid.nn<-pred.resid.nn*max2.resp.f
				pred.resid.nn<-pred.resid.nn+min.resp.f
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-dat_tps[[i]]$resp-pred.resid.nn
				res.FINAL<-res.FINAL+(res.FINAL.2*OptX.mfit.wt[[iter.mod]])}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(dat_tps[[i]]$resp-pred.resid.nn)*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
             			
  	          ##final models = boosted regression tree	
			  if (k=="b"){
				gc()
			  	print(paste("Final modeling of ",(out.names[i]),": Boosted Regression Trees"))
				print(paste("                                         ",Sys.time()))
				mod.run="BRT"
				#run model with all points
				mod.brt.tps.FINAL<- machisplin.gbm.step(data=dat_tps[[i]], gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5, plot.main = FALSE, silent = TRUE)
				#store variable importance
				l$var.imp$brt<-mod.brt.tps.FINAL$contributions
				#create raster prediction
				if(is.null(pred.elev)==FALSE){pred.elev.2<- terra::predict(rast_stack, mod.brt.tps.FINAL, n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")
				pred.elev<-pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]])}
				if(is.null(pred.elev)==TRUE){pred.elev<-terra::predict(rast_stack, mod.brt.tps.FINAL, n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")*OptX.mfit.wt[[iter.mod]]}
				#predict at all train sites
				pred.brt.obs <- terra::predict(mod.brt.tps.FINAL, dat_tps[[i]], n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")
				#get residuals
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-(dat_tps[[i]]$resp-pred.brt.obs)*OptX.mfit.wt[[iter.mod]]
				res.FINAL<-res.FINAL+res.FINAL.2}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(dat_tps[[i]]$resp-pred.brt.obs)*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
			  
  		      ##final models = randomForest	
         	  #if (mfit.rf<mfit.brt & mfit.rf<mfit.lm & mfit.rf=<mfit.mars){
			  if (k=="r"){
				gc()
				mod.run="RF"
				print(paste("Final modeling of ",(out.names[i]),": Random Forests"))
				print(paste("                                         ",Sys.time()))
				#run model with all points
				mod.rf.tps.FINAL<- randomForest::randomForest(mod.form, data = dat_tps[[i]],importance = TRUE)
				#store variable importance
				l$var.imp$rf<-mod.rf.tps.FINAL$importance
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-terra::predict(rast_stack, mod.rf.tps.FINAL, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
				pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(pred.elev)==TRUE){pred.elev<-terra::predict(rast_stack, mod.rf.tps.FINAL, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				pred.rf.obs<-terra::predict(mod.rf.tps.FINAL, dat_tps[[i]], type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-(dat_tps[[i]]$resp-pred.rf.obs)*OptX.mfit.wt[[iter.mod]]
				res.FINAL<-(res.FINAL+res.FINAL.2)}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(dat_tps[[i]]$resp-pred.rf.obs)*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
                
			  ##final models = MARS	
			  if (k=="m"){
				gc()
				mod.run="MARS"
				print(paste("Final modeling of ",(out.names[i]),": Multivariate Adaptive Regression Splines"))
				print(paste("                                         ",Sys.time()))
				#run model with all points
				mod.MARS.tps.FINAL<-earth::earth(mod.form,  data = dat_tps[[i]], nfold=10)
				#store variable importance
				l$var.imp$mars<-earth::evimp(mod.MARS.tps.FINAL)
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-terra::predict(rast_stack, mod.MARS.tps.FINAL)
				pred.elev<-pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]])}
				if(is.null(pred.elev)==TRUE){pred.elev<-terra::predict(rast_stack, mod.MARS.tps.FINAL)*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-resid(mod.MARS.tps.FINAL, warn=FALSE)
				res.FINAL<-(res.FINAL+(as.vector(res.FINAL.2)*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(as.vector(resid(mod.MARS.tps.FINAL, warn=FALSE)))*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
              
   		      ##final models = SVM
			  if (k=="v"){
				gc()
				mod.run="SVM"
				print(paste("Final modeling of ",(out.names[i]),": Support Vector Machines"))
				print(paste("                                         ",Sys.time()))
				#run model with all points
				mod.SVM.tps.FINAL<-kernlab::ksvm(mod.form, data=dat_tps[[i]])
				#store variable importance
				nvals.in<-nrow(dat_tps[[i]])
				sampleSVM<-dat_tps[[i]]
				if(nrow(sampleSVM)>200){sampleSVM<- sampleSVM[sample(nrow(sampleSVM),200),]}
				nvals.in<-nrow(sampleSVM)
				SVM.imp<-NULL
				for(sub.svm in 1:nvals.in){
				nob <- sampleSVM[sub.svm,]
				nob <- nob[1:(1+n.covars)]
				set.seed(1313)	
				explain_z<- breakDown::broken(mod.SVM.tps.FINAL,  new_observation = nob, data = sampleSVM, predict.function = terra::predict, baseline = "intercept", direction = "up", ) 
				if(is.null(SVM.imp)==TRUE){SVM.imp<-abs(explain_z$contribution)}
				if(is.null(SVM.imp)==FALSE){SVM.imp<-abs(SVM.imp)+abs(explain_z$contribution)}
				}
				SVM.f<-(SVM.imp/nvals.in)
				SVM.f<-as.matrix(SVM.f,ncol=1)
				name.SVM.imp<-gsub(" =.*","",explain_z$variable)
				rownames(SVM.f)<-name.SVM.imp
				colnames(SVM.f)<-("contributions to SVM")
				l$var.imp$SVM<-SVM.f
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-terra::predict(rast_stack, mod.SVM.tps.FINAL, na.rm=TRUE)
				pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(pred.elev)==TRUE){pred.elev<-terra::predict(rast_stack, mod.SVM.tps.FINAL, na.rm=TRUE)*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				pred.SVM.obs<-terra::predict(mod.SVM.tps.FINAL, dat_tps[[i]])
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-as.vector(dat_tps[[i]]$resp-pred.SVM.obs)
				res.FINAL<-(res.FINAL+(res.FINAL.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(as.vector(dat_tps[[i]]$resp-pred.SVM.obs))*OptX.mfit.wt[[iter.mod]]}
				gc()
				}
              
   		      ##final models = GAM
			  if (k=="g"){
				gc()
				mod.run="GAM"
				print(paste("Final modeling of ",(out.names[i]),": Generalized Additive Model"))				
				print(paste("                                         ",Sys.time()))
				#run model with all points
				mod.GAM.tps.FINAL<-mgcv::gam(mod.form, data=dat_tps[[i]])
				#store variable importance
				l$var.imp$gam<-mod.GAM.tps.FINAL$coefficients
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-terra::predict(rast_stack, mod.GAM.tps.FINAL)
				pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(pred.elev)==TRUE){pred.elev<-(terra::predict(rast_stack, mod.GAM.tps.FINAL))*OptX.mfit.wt[[iter.mod]]}
				#get residuals
				pred.GAM.obs<-terra::predict(mod.GAM.tps.FINAL, dat_tps[[i]])
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-as.vector(dat_tps[[i]]$resp-pred.GAM.obs)
				res.FINAL<-(res.FINAL+(res.FINAL.2*OptX.mfit.wt[[iter.mod]]))}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-(as.vector(dat_tps[[i]]$resp-pred.GAM.obs))*OptX.mfit.wt[[iter.mod]]}
	   		    gc()
				}
			}
			print(paste("Calculating final ensemble of models for ",(out.names[i])))
			print(paste("                                         ",Sys.time()))
			#wrap up analysis and get final layer
			#divide models by number ensembled
            pred.elev<-(pred.elev/OptX.mfit.wt.tot)
			res.FINAL<-(res.FINAL/OptX.mfit.wt.tot)
			
			#calculate Sum of squares and resisdual sum of squares
			print(paste("Calculating residuals of final model ensemble: ",(out.names[i])))
			print(paste("                                         ",Sys.time()))
			rss.m <- sum((res.FINAL) ^ 2)
			tss <- sum((dat_tps[[i]]$resp - mean(dat_tps[[i]]$resp)) ^ 2)
			l$residuals<- cbind((res.FINAL),dat_tps[[i]]$LONG,dat_tps[[i]]$LAT)
			colnames(l$residuals)<-c("residuals","long","lat")

			rsq.model <- 1 - (rss.m/tss) #r squared
			gc()
			##################################################################################################
			###################################### part 3 spline residuals ###################################
			##################################################################################################
			#DETERMINE IF TILING IS NEEDED
			if(tps==TRUE) {
				print("###########################################################################################")
				print("###########################################################################################")
				print("##################################### STEP 3 (of 3) #######################################")
				print("###########################################################################################")
				print("######### Correcting model error: thin plate spline of residuals of #######################")
				print("###########################################################################################")
				print(paste("###########################  ", (out.names[i]), "  ######################################"))
				print("###########################################################################################")
				print("###########################################################################################")
				print(paste("                                         ",Sys.time()))
				options(warn = -1)
				#specify total extents
				totalExt<-terra::ext(rast_stack)
				#specify n rows 
				nr.in<-nrow(rast_stack)
				#specify n cols 
				nc.in<-ncol(rast_stack)
				#specify n row blocks,3000
				#nrB<-nr.in/3000
				nrB<-nr.in/3000
				nRx<-ceiling(nrB)
				#specify n col blocks,3000
				#ncB<-nc.in/3000
				ncB<-nc.in/3000
				nCx<-ceiling(ncB)
				#specify lat distance
				longDist<-((totalExt[2]-totalExt[1])/nCx)
				#specify long distance
				latDist<-((totalExt[4]-totalExt[3])/nRx)
				
				#specify sampling grid for TPS
				m2<-0;m1<-0;new.df<-NULL;new.df2<-NULL
				
				for (j in 1:nRx){
				for (h in 1:nCx){
						m1<-m1+1
						new.df[[m1]]<-c((totalExt[1]+((longDist*(h-1))-(longDist*0.2))),(totalExt[1]+((longDist*h)+(longDist*0.2))),(totalExt[3]+((latDist*(j-1)))-(latDist*0.2)),(totalExt[3]+((latDist*j))+(latDist*0.2)))
					}
				}
				#specify sampling grid for mosaic
				for (j in 1:nRx){
				   for (h in 1:nCx){
						m2<-m2+1
						new.df2[[m2]]<-c((totalExt[1]+((longDist*(h-1))-(longDist*0.025))),(totalExt[1]+((longDist*h)+(longDist*0.025))),(totalExt[3]+((latDist*(j-1)))-(latDist*0.025)),(totalExt[3]+((latDist*j))+(latDist*0.025)))
				}}
				#specify full cordinates for spatial subsampling
				#perform TPS on each grid and then mosaic
				print(paste("Thin plate splines of residuals will be tiled across ",(nRx*nCx)," tile(s), GIS layer: ", (out.names[i])))
				print(paste("                                         ",Sys.time()))
				#perform TPS on each grid and then mosaic
				if(nRx*nCx>1){
					Full.cords<-dat_tps[[1]][,c(n.covars,n.covars+1)]
					pred_TPS_elev<-NULL
					for (h in 1:(nRx*nCx)){
						gc()
						print(paste("Performing thin plate splines of residuals on tile",(h)))
						print(paste("                                         ",Sys.time()))
						#clip raster brick, CURRENT WORKING DIRECTORY 9/10
						b <- terra::ext(new.df[[h]][1], new.df[[h]][2], new.df[[h]][3], new.df[[h]][4])
						d <- terra::ext(new.df2[[h]][1], new.df2[[h]][2], new.df2[[h]][3], new.df2[[h]][4])
						#terra::crs(b) <- terra::crs(rast_stack)
						#terra::crs(d) <- terra::crs(rast_stack)
						rb <- terra::crop(rast_stack, b)
						#sub sample residuals
						RAST_TPS<-data.frame(terra::extract(rb[[1]], Full.cords))
						#str(RAST_VAL)
						#merge sampled data to input
						MyTPSdata<-cbind(res.FINAL,Full.cords,RAST_TPS)
						#remove NAs
						MyTPSdata<-MyTPSdata[complete.cases(MyTPSdata),][1:3]
						
						print(paste("Number of data points in grid:",nrow(MyTPSdata)))
						#fit thin plate spline of residuals
						if(nrow(MyTPSdata)<10){
							print("Sample grid mostly empty: grid area likely is water or a region with no input data. This grid will not be part of TPS. This is expected in many situations and part of the error correction process.")
							rbT<-rb[[1]]
							r<-terra::rast(rbT)
							r[is.na(r[])] <- 0
							pred_TPS_elev<-r
							TPS_name<-paste0("TILE_",h)
							names(pred_TPS_elev)<- TPS_name
							pred_TPS_elev <- terra::crop(pred_TPS_elev, d)
							pred_TPS_elev <- terra::extend(pred_TPS_elev,rast_stack)
							terra::plot(pred_TPS_elev)
							Sys.sleep(0)} 
						else {mod.tps.elev<-fields::Tps(MyTPSdata[2:3], MyTPSdata[1])#columns of lat and long
							#use TPS to interpolate residuals
							TPS_name<-paste0("TILE_",h)
							#use TPS to interpolate residuals
							pred_TPS_elev<-terra::interpolate(terra::rast(rb), mod.tps.elev)
							names(pred_TPS_elev)<- TPS_name
							pred_TPS_elev <- terra::crop(pred_TPS_elev, d)
							pred_TPS_elev <- terra::extend(pred_TPS_elev,rast_stack)
							terra::plot(pred_TPS_elev)
							Sys.sleep(0)}
						if (h>1){
							r_list <- list(raster_sTPS, pred_TPS_elev)
							raster_sTPS <- terra::rast(r_list)
							pred_TPS_elev<-NULL}
						else {raster_sTPS<-pred_TPS_elev}
						gc()
						}}
				if(nRx*nCx>1){
						for (j in 1:terra::nlyr(raster_sTPS)) {
							    if(j==1){
								   ie<-terra::sprc(raster_sTPS[[1]])
								}else{
								   ie<-terra::sprc(raster_sTPS[[j]],ie)
								}}
						rast.mosaic<- terra::mosaic(ie, fun= "mean")
								}
				if(nRx*nCx==1){
					#fit thin plate spline of residuals
					print("Performing thin plate splines of residuals")				
					mod.tps.elev<-fields::Tps(dat_tps[[i]][,c(n.covars,n.covars+1)], res.FINAL)#columns of lat and long
					#use TPS to interpolate residuals
					rast.mosaic<-terra::interpolate(terra::rast(rast_stack), mod.tps.elev)}
				
				##################################################################################################
				############################## part 4 feather seams of spline tiles ##############################
				##################################################################################################
				#feather right
				feath.ras.TPS<-NULL
				if(nRx*nCx>1){
					print("Feathering seams of thin plate splines")
					print(paste("Feathering",nRx,"rows and",nCx,"columns"))
					print(paste("                                         ",Sys.time()))
					for (j in 1:nRx){
					  for (h in 1:nCx){
						gc()
						v<-(h+((j*nCx)-nCx))
						if (h<nCx){
						print(paste("Feathering vertical seams of thin-plate splines for row",j,"and column",h))		
						print(paste("Seam is at edge of",v,"and",v+1,"tiles"))
						print(paste("                                         ",Sys.time()))
						AAA<-raster_sTPS[[v]]+raster_sTPS[[v+1]]
						MyShp<-terra::as.points(AAA, value=TRUE)						
						#specify coordinate system
						terra::crs(MyShp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
						#
						blend.ext<-terra::ext(MyShp)
						#crop
						bR1 <- terra::crop(raster_sTPS[[v]], blend.ext)
						bR2 <- terra::crop(raster_sTPS[[v+1]], blend.ext)
						##
						x_bR1<- terra::xFromCell(bR1, 1:terra::ncell(bR1))
						x_bR2<- terra::xFromCell(bR2, 1:terra::ncell(bR2))
						
						#creating fade
						delta.L1<-(max(x_bR1)-(min(x_bR1)))
						stD1<-(x_bR1-min(x_bR1))/delta.L1
						stD1<-1-stD1
						delta.L2<-(max(x_bR2)-(min(x_bR2)))
						stD2<-(x_bR1-min(x_bR2))/delta.L2
						#latitude_values <- terra::yFromCell(r, 1:terra::ncell(r))LONG.B1 <- bR1
						
						we.R.d1a <- terra::setValues(bR1, stD1)
						we.R.d2a <- terra::setValues(bR2, stD2)
						we.R.d1 <-bR1 * we.R.d1a
						we.R.d2 <-bR2 * we.R.d2a
						feath.ras<-we.R.d2+we.R.d1
						feath.ras<- terra::extend(feath.ras,rast_stack)
						terra::plot(feath.ras)
						if(is.null(feath.ras.TPS)==TRUE){
							feath.ras.TPS<-feath.ras
						}else {
						r_list <- list(feath.ras, feath.ras.TPS)
						feath.ras.TPS<- terra::rast(r_list)
						}}}}
						gc()
						#####################################################
						#####################################################
					print("Finished feathering vertical edges")
					print("##################################")
					print("Starting feathering horizontal edges")
					print(paste("Feathering",nRx,"rows and",nCx,"columns"))
					f.clock=0
					for (j in 1:nRx){
					  for (h in 1:nCx){
						#feather up
						f.clock=f.clock+1
						f.timer=(nRx*nCx)-nCx+1
						#v<-(h+((j*nRx)-nRx))
						v<-(h+((j*nCx))-nCx)
						g<-(h+((j*nCx)))
						if (f.clock<f.timer){
						print(paste("Feathering horizonal seams of thin-plate splines for row ",j,"and column",h))	
						print(paste("Seam is at edge of",v,"and",g,"tiles"))
						print(paste("                                         ",Sys.time()))
						gc()
						raster_sTPS<-raster_sTPS
						AAA<-raster_sTPS[[v]]+raster_sTPS[[g]]
						MyShp<-terra::as.points(AAA, value=TRUE)
						terra::crs(MyShp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
						#
						blend.ext<-terra::ext(MyShp)
						#crop
						bR1 <- terra::crop(raster_sTPS[[v]], blend.ext)
						bR2 <- terra::crop(raster_sTPS[[v+nCx]], blend.ext)
						if(identical(bR1,bR2)==FALSE){
							#trim.ext<-terra::ext(we.R.d1a)
							bR1.Rin<-nrow(bR1)
							bR1.Cin<-ncol(bR1)
							bR2.Rin<-nrow(bR2)
							bR2.Cin<-ncol(bR2)
							if(bR1.Rin>bR2.Rin){
								fix.ext<-terra::ext(bR2)
								bR1 <- terra::crop(bR1, fix.ext)
								bR1 <- terra::crop(bR1, fix.ext)
							  }else{
							    fix.ext<-terra::ext(bR1)
								bR2 <- terra::crop(bR2, fix.ext)
							   }
							}						
						##sample long and lat
						x_bR1<- terra::yFromCell(bR1, 1:terra::ncell(bR1))
						x_bR2<- terra::yFromCell(bR2, 1:terra::ncell(bR2))
						#creating fade
						delta.L1<-(max(x_bR1)-(min(x_bR1)))
						stD1<-(x_bR1-min(x_bR1))/delta.L1
						stD1<-1-stD1
						delta.L2<-(max(x_bR2)-(min(x_bR2)))
						stD2<-(x_bR1-min(x_bR2))/delta.L2
						
						we.R.d1a <- terra::setValues(bR1, stD1)
						we.R.d2a <- terra::setValues(bR2, stD2)
						we.R.d1 <-bR1 * we.R.d1a
						we.R.d2 <-bR2 * we.R.d2a
						
						feath.ras<-we.R.d2+we.R.d1
						feath.ras<- terra::extend(feath.ras,rast_stack)
						terra::plot(feath.ras)
						if(is.null(feath.ras.TPS)==TRUE){
							feath.ras.TPS<-feath.ras
						} else {
						r_list <- list(feath.ras, feath.ras.TPS)
							feath.ras.TPS<- terra::rast(r_list)
							}
						gc()
						}}}
					print("Finished feathering horizontal edges")
					print("Merging all raster tiles into single raster")
					if(nRx*nCx>2){
						for (j in 1:terra::nlyr(feath.ras.TPS)) {
							    if(j==1){
								   ic<-terra::sprc(feath.ras.TPS[[1]])
								}else{
								   ic<-terra::sprc(feath.ras.TPS[[j]],ic)
								}}
						rast.mosaic.TPS <- terra::mosaic(ic, fun= "mean")
						etal<-terra::sprc(rast.mosaic.TPS,rast.mosaic)
						final.TPS<- terra::mosaic(etal, fun= "first")
						}
					if(nRx*nCx==2){
						final.TPS<-terra::merge(feath.ras.TPS,rast.mosaic)}
					if(nRx*nCx==1){
						final.TPS<-rast.mosaic
						}
			}
			}
			print("Finalizing Results")
			##################################################################################################
			################################# part 5 finalize results ##################################
			##################################################################################################
			if(tps==TRUE) {
				gc()
				print("Creating final data files")
				#calculate pred at normal scale, suming up kriging pred+res
				pred.elev.i<-c(pred.elev, final.TPS)
				pred.elev.i.calc<- terra::app(pred.elev.i, fun=sum)
				
				#extract final values to input points from final raster
				f.actual<-terra::extract(pred.elev.i.calc, dat_tps[[i]][,n.covars:(n.covars+1)])#lat and long input
            
				#calculate Sum of squares and resisdual sum of squares
				rss.final <- sum((dat_tps[[i]]$resp - f.actual[2]) ^ 2)
				l$residuals<- cbind((dat_tps[[i]]$resp - f.actual[2]),dat_tps[[i]]$LONG,dat_tps[[i]]$LAT)
				colnames(l$residuals)<-c("residuals","long","lat")
				#tss <- sum((dat_tps[[i]]$resp - mean(dat_tps[[i]]$resp)) ^ 2)
				rsq.final <- 1 - (rss.final/tss) #r squared
				
				#out values
				sum.val<-data.frame(out.names[i],f.mod,OptX.mfit.txt,rsq.model,rsq.final)
				colnames(sum.val)<-c("layer","best model(s):","ensemble weights:","r2 ensemble:","r2 final:")
				l$n.layers<-length(out.names)
				l$summary<-sum.val
				#save TPS if R2 is better, else save only model
				if(rsq.final>rsq.model){
					names(pred.elev.i.calc)<- out.names[i]
					l$final<-pred.elev.i.calc} else {
					names(pred.elev)<- out.names[i]
					l$final<-pred.elev
					}
			options(warn = 0)
			}
			
			if(tps==FALSE) {
				print("###########################################################################################")
				print("###########################################################################################")
				print("##################################### STEP 3 (of 3) #######################################")
				print("###########################################################################################")
				print("######### SKIPPING: Correcting model error via thin-plate splines of residuals #############")
				print("###################  You selected input parameter: tps= FALSE #############################")
				print(paste("###########################  ", (out.names[i]), "  ######################################"))
				print("###########################################################################################")
				print("###########################################################################################")
				gc()
				#out values
				sum.val<-data.frame(out.names[i],f.mod,OptX.mfit.txt,rsq.model)
				colnames(sum.val)<-c("layer","best model(s):","ensemble weights:","r2 ensemble:")
				if(i==1){l$n.layers<-length(out.names)}
				l$summary<-sum.val
				#save TPS if R2 is better, else save only model
				names(pred.elev)<- out.names[i]
				l$final<-pred.elev
				}
		#l$var.imp
		omega[[i]] <-l
		print(paste("Finalizing All Results!!!:", (out.names[i]), Sys.time()))
		print("*******************************************************************************************")
		print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
		print("###########################################################################################")
		print(paste("############################ FINSHED:",(out.names[i]),"###############################"))
		print("###########################################################################################")
		print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
		print("*******************************************************************************************")
		
	}
sink()
return(omega)}
}

##################################################################################################
##################################################################################################
#' Write MACHISPLIN geotiff
#' @param mttps.in an output from the 'machisplin.mltps' function  
#' @param out.names a vector corresponding to the output raster names.   If 'null' it will write intial input names.  
#' @return This function outputs performance values, the algorithm(s) used, and rasters for use in GIS output from the 'machisplin.mltps' function
#' @export
#' @import terra
#' @examples
#' library(MACHISPLIN)
#' 
#' # Import a csv 
#' Mydata	<-read.delim("sampling.csv", sep=",", h=T)
#' 
#' #load rasters to use as high resolution co-variates for downscaling
#' ALT = raster("SRTM30m.tif")
#' SLOPE = raster("ln_slope.tif")
#' ASPECT = raster("aspect.tif")
#' TWI = raster("TWI.tif")
#' 
#' # function input: raster stack of covarites
#' raster_stack<-stack(ALT,SLOPE,TWI,ASPECT)
#' 
#' #run an ensemble machine learning thin plate spline 
#' interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2)
#' 
#' machisplin.write.geotiff(mltps.in=interp.rast)

machisplin.write.geotiff<-function(mltps.in, out.names= NULL, overwrite=TRUE){
	n.spln<-mltps.in[[1]]$n.layers
	
	### store rasters
    for (i in 1:n.spln){
		if(i ==1){rast.O<-mltps.in[[i]]$final} else {rast.O<-c(rast.O, mltps.in[[i]]$final)}
		}
	output_dir <- getwd()
	#write rasters
	if(is.null(out.names)==TRUE){
	   for (i in 1:n.spln) {
		layer_name <-names(rast.O[[i]])
		output_filename <- file.path(output_dir, paste0(layer_name, ".tif"))
		terra::writeRaster(rast.O[[i]], filename = output_filename, overwrite = TRUE)
		print(paste("Layer", layer_name, "saved to:", output_filename))
		}
	}
	
	if(is.null(out.names)==FALSE){
	   for (i in 1:n.spln) {
		layer_name <-out.names[i]
		output_filename <- file.path(output_dir, paste0(layer_name, ".tif"))
		terra::writeRaster(rast.O[[i]], filename = output_filename, overwrite = TRUE)
		print(paste("Layer", layer_name, "saved to:", output_filename))
		}
	}
	#if(is.null(out.names)==TRUE){terra::writeRaster(c(rast.O), names(rast.O), format="GTiff",overwrite=overwrite)}
	#if(is.null(out.names)==FALSE){terra::writeRaster(c(rast.O), names(out.names), format="GTiff",overwrite=overwrite)}

	## make summary table - add time and date to output as a row
	for (i in 1:n.spln){
		if(i ==1){table.O<-mltps.in[[i]]$summary} else{table.O<-rbind(table.O,mltps.in[[i]]$summary)}
		}
	#write summary table
	write.csv(table.O, "MACHISPLIN_results.csv")
	
	legend<-""
	legend0<-"R2 Final: ensemble of the best models & thin-plate-spline of the residuals of the ensemble model"
	legend1<-"Best model legend: The quantity of letters depicts the number of models ensembled." 
	legend2<-"The letters themselves depict the model algorithm: b = boosted regression trees (BRT);"
	legend3<-"g = generalized additive model (GAM); m = multivariate adaptive regression splines (MARS);"
	legend4<-"v = support vector machines (SVM); r = random forests (RF); n = neural networks (NN)"
	legend5<-"The ensemble weights is percentage that each algorithm contributed to the ensemble model"
	legend6<-"NOTE: if 'R2 Ensemble' is greater than 'R2 Final', then the output model is only the ensembled model (the thin-plate-spline of residuals were not used)"
	write(legend,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend0,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend1,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend2,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend3,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend4,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend5,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend6,file="MACHISPLIN_results.csv",append=TRUE)
}

##################################################################################################
##################################################################################################
#' Write MACHISPLIN model loadings
#' @param mttps.in an output from the 'machisplin.mltps' function  
#' @param out.names a vector corresponding to the output raster names.   If 'null' it will write intial input names.  
#' @return This function outputs performance values, the algorithm(s) used, and rasters for use in GIS output from the 'machisplin.mltps' function
#' @export
#' @import terra
#' @examples
#' library(MACHISPLIN)
#' 
#' # Import a csv 
#' Mydata	<-read.delim("sampling.csv", sep=",", h=T)
#' 
#' #load rasters to use as high resolution co-variates for downscaling
#' ALT = raster("SRTM30m.tif")
#' SLOPE = raster("ln_slope.tif")
#' ASPECT = raster("aspect.tif")
#' TWI = raster("TWI.tif")
#' 
#' # function input: raster stack of covarites
#' raster_stack<-stack(ALT,SLOPE,TWI,ASPECT)
#' 
#' #run an ensemble machine learning thin plate spline 
#' interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2)
#' 
#' machisplin.write.loadings(mltps.in=interp.rast)
machisplin.write.loadings<-function(mltps.in, out.names= NULL){
	n.spln<-mltps.in[[1]]$n.layers
	for (i in 1:n.spln){
		sink(paste0(mltps.in[[i]]$summary[[1]],"_model_loadings.txt"))
		print(mltps.in[[i]]$var.imp)	
		}
sink()	
}
 # returns output to the console
##################################################################################################
##################################################################################################
#' Write MACHISPLIN model residuals
#' @param mttps.in an output from the 'machisplin.mltps' function  
#' @param out.names a vector corresponding to the output raster names.   If 'null' it will write intial input names.  
#' @return This function outputs performance values, the algorithm(s) used, and rasters for use in GIS output from the 'machisplin.mltps' function
#' @export
#' @import terra
#' @examples
#' library(MACHISPLIN)
#' 
#' # Import a csv 
#' Mydata	<-read.delim("sampling.csv", sep=",", h=T)
#' 
#' #load rasters to use as high resolution co-variates for downscaling
#' ALT = raster("SRTM30m.tif")
#' SLOPE = raster("ln_slope.tif")
#' ASPECT = raster("aspect.tif")
#' TWI = raster("TWI.tif")
#' 
#' # function input: raster stack of covarites
#' raster_stack<-stack(ALT,SLOPE,TWI,ASPECT)
#' 
#' #run an ensemble machine learning thin plate spline 
#' interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2)
#' 
#' machisplin.write.residuals(mltps.in=interp.rast)

machisplin.write.residuals<-function(mltps.in, out.names= NULL){
	n.spln<-mltps.in[[1]]$n.layers
	for (i in 1:n.spln){
        table.O<-mltps.in[[i]]$residuals
		write.csv(table.O, paste0(mltps.in[[i]]$summary[[1]],"_residuals.csv"))
		}
}

##################################################################################################
##################################################################################################
#library(terra);library(gstat);library(MASS);library(nlme);library(fields);library(randomForest);library(earth);library(kernlab);library(mgcv);library(snowfall);library(snow);library(spdep);library(nnet);library(optimx);library(NeuralNetTools);library(breakDown)



machisplin.kfold <- function(x, k=5, by=NULL) {

	singlefold <- function(obs, k) {
		if (k==1) {
			return(rep(1, obs))
		} else {
			i <- obs / k
			if (i < 1) {
				stop('insufficient records:', obs, ', with k=', k)
			}
			i <- round(c(0, i * 1:(k-1), obs))
			times = i[-1] - i[-length(i)]

			group <- c()
			for (j in 1:(length(times))) {
				group <- c( group, rep(j, times=times[j]) )
			}
		
			r <- order(runif(obs))
			return(group[r]) 
		}
	}

	if (is.vector(x)) {
		if (length(x) == 1) {
			if (x > 1) {
				x <- 1:x
			}
		}
		obs <- length(x)
	} else if (inherits(x, 'Spatial')) {
		if (inherits(x, 'SpatialPoints')) {
			obs <- nrow(coordinates(x))
		} else {
			obs <- nrow(x@data)
		}
	} else {
		obs <- nrow(x)
	}
	if (is.null(by)) {
		return(singlefold(obs, k))
	}
	
	by = as.vector(as.matrix(by))
	if (length(by) != obs) {
		stop('by should be a vector with the same number of records as x')
	}
	un <- unique(by)
	group <- vector(length=obs)
	for ( u in un ) {
		i = which(by==u)
		kk = min(length(i), k)
		if (kk < k) warning('lowered k for by group: ', u  ,'  because the number of observations was  ',  length(i))
		group[i] <- singlefold(length(i), kk)
	} 
	return(group)
}


#
# j. leathwick/j. elith - 19th September 2005
#
# version 2.9
#
# function to assess optimal no of boosting trees using k-fold cross validation 
#
# implements the cross-validation procedure described on page 215 of 
# Hastie T, Tibshirani R, Friedman JH (2001) The Elements of Statistical Learning: 
# Data Mining, Inference, and Prediction Springer-Verlag, New York.
#
# divides the data into 10 subsets, with stratification by prevalence if required for pa data
# then fits a gbm model of increasing complexity along the sequence from n.trees to n.trees + (n.steps * step.size)
# calculating the residual deviance at each step along the way
# after each fold processed, calculates the average holdout residual deviance and its standard error
# then identifies the optimal number of trees as that at which the holdout deviance is minimised 
# and fits a model with this number of trees, returning it as a gbm model along with additional information 
# from the cv selection process
#
# updated 13/6/05 to accommodate weighting of sites
#
# updated 19/8/05 to increment all folds simultaneously, allowing the stopping rule
# for the maxinum number of trees to be fitted to be imposed by the data, 
# rather than being fixed in advance
#
# updated 29/8/05 to return cv test statistics, and deviance as mean
# time for analysis also returned via unclass(Sys.time())
#
# updated 5/9/05 to use external function calc.deviance 
# and to return cv test stats via predictions formed from fold models 
# with n.trees = target.trees
#
# updated 15/5/06 to calculate variance of fitted and predicted values across folds
# these can be expected to approximate the variance of fitted values
# as would be estimated for example by bootstrapping
# as these will underestimate the true variance 
# they are corrected by multiplying by (n-1)2/n 
# where n is the number of folds
#
# updated 25/3/07 tp allow varying of bag fraction
#
# requires gbm library from Cran
# requires roc and calibration scripts of J Elith
# requires calc.deviance script of J Elith/J Leathwick
#
#


machisplin.gbm.step <- function (
  data,                                     # the input dataframe
  gbm.x,                                    # the predictors
  gbm.y,                                    # and response
  offset = NULL,                            # allows an offset to be specified
  fold.vector = NULL,                       # allows a fold vector to be read in for CV with offsets,
  tree.complexity = 1,                      # sets the complexity of individual trees
  learning.rate = 0.01,                     # sets the weight applied to inidivudal trees
  bag.fraction = 0.75,                      # sets the proportion of observations used in selecting variables
  site.weights = rep(1, nrow(data)),        # allows varying weighting for sites
  var.monotone = rep(0, length(gbm.x)),     # restricts responses to individual predictors to monotone 
  n.folds = 10,                             # number of folds
  prev.stratify = TRUE,                     # prevalence stratify the folds - only for p/a data
  family = "bernoulli",                     # family - bernoulli (=binomial), poisson, laplace or gaussian
  n.trees = 50,                             # number of initial trees to fit
  step.size = n.trees,                      # numbers of trees to add at each cycle
  max.trees = 10000,                        # max number of trees to fit before stopping
  tolerance.method = "auto",                # method to use in deciding to stop - "fixed" or "auto"
  tolerance = 0.001,                        # tolerance value to use - if method == fixed is absolute, 
                                            # if auto is multiplier * total mean deviance
  plot.main = TRUE,                         # plot hold-out deviance curve
  plot.folds = FALSE,                       # plot the individual folds as well
  verbose = TRUE,                           # control amount of screen reporting
  silent = FALSE,                           # to allow running with no output for simplifying model)
  keep.fold.models = FALSE,                 # keep the fold models from cross valiation
  keep.fold.vector = FALSE,                 # allows the vector defining fold membership to be kept
  keep.fold.fit = FALSE,                    # allows the predicted values for observations from CV to be kept
  ...)                                      # allows for any additional plotting parameters
{

	if (! requireNamespace('gbm') ) { stop ('you need to install the gbm package to run this function') }
	requireNamespace('splines')

	if (silent) verbose <- FALSE

# initiate timing call

	z1 <- Sys.time()

# setup input data and assign to position one

#  dataframe.name <- deparse(substitute(data))   # get the dataframe name
#  data <- eval(data)
	x.data <- data[, gbm.x, drop=FALSE]                 #form the temporary datasets
	#names(x.data) <- names(data)[gbm.x]
	y.data <- data[, gbm.y]
	sp.name <- names(data)[gbm.y]
	if (family == "bernoulli") {
		prevalence <- mean(y.data)
	}

#  assign("x.data", x.data, env = globalenv())               #and assign them for later use
#  assign("y.data", y.data, env = globalenv())

	#offset.name <- deparse(substitute(offset))   # get the dataframe name
	#offset <- eval(offset)

	n.cases <- nrow(data)
	n.preds <- length(gbm.x)

	if (!silent) {
		cat("\n","\n","GBM STEP - version 2.9","\n","\n")
		cat("Performing cross-validation optimisation of a boosted regression tree model \n")
		cat("for", sp.name, "and using a family of",family,"\n")
		cat("Using",n.cases,"observations and",n.preds,"predictors \n")
	}

# set up the selector variable either with or without prevalence stratification

	if (is.null(fold.vector)) {

		if (prev.stratify & family == "bernoulli") {
			presence.mask <- data[,gbm.y] == 1
			absence.mask <- data[,gbm.y] == 0
			n.pres <- sum(presence.mask)
			n.abs <- sum(absence.mask)

# create a vector of randomised numbers and feed into presences
			selector <- rep(0, n.cases)
			temp <- rep(seq(1, n.folds, by = 1), length = n.pres)
			temp <- temp[order(runif(n.pres, 1, 100))]
			selector[presence.mask] <- temp

# and then do the same for absences
			temp <- rep(seq(1, n.folds, by = 1), length = n.abs)
			temp <- temp[order(runif(n.abs, 1, 100))]
			selector[absence.mask] <- temp
		} else {  #otherwise make them random with respect to presence/absence
			selector <- rep(seq(1, n.folds, by = 1), length = n.cases)
			selector <- selector[order(runif(n.cases, 1, 100))]
		}
    } else {
		if (length(fold.vector) != n.cases) {
			stop("supplied fold vector is of wrong length")
		}
		cat("loading user-supplied fold vector \n")
		selector <- fold.vector
    }

# set up the storage space for results

	pred.values <- rep(0, n.cases)

	cv.loss.matrix <- matrix(0, nrow = n.folds, ncol = 1)
	training.loss.matrix <- matrix(0, nrow = n.folds, ncol = 1)
	trees.fitted <- n.trees

	model.list <- list(paste("model", 1:n.folds, sep=""))     # dummy list for the tree models

# set up the initial call to gbm

	if (is.null(offset)) {
		gbm.call <- paste("gbm::gbm(y.subset ~ .,data=x.subset, n.trees = n.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = weight.subset, distribution = as.character(family), var.monotone = var.monotone, verbose = FALSE)", sep="")
    } else {
		gbm.call <- paste("gbm::gbm(y.subset ~ . + offset(offset.subset), data=x.subset, n.trees = n.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = weight.subset, distribution = as.character(family), var.monotone = var.monotone, verbose = FALSE)", sep="")
    } 

	n.fitted <- n.trees

# calculate the total deviance
               
	y_i <- y.data

	u_i <- sum(y.data * site.weights) / sum(site.weights)
	u_i <- rep(u_i,length(y_i))

	total.deviance <- machisplin.calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)

	mean.total.deviance <- total.deviance/n.cases

	tolerance.test <- tolerance

	if (tolerance.method == "auto") {
		tolerance.test <- mean.total.deviance * tolerance
	}

# now step through the folds setting up the initial call

	if (!silent){ 
		cat("creating",n.folds,"initial models of",n.trees,"trees","\n")
		if (prev.stratify & family == "bernoulli") {
			cat("\n","folds are stratified by prevalence","\n")
		} else {
			cat("\n","folds are unstratified","\n")
		}
		cat ("total mean deviance = ",round(mean.total.deviance,4),"\n")
		cat("tolerance is fixed at ",round(tolerance.test,4),"\n")
		if (tolerance.method != "fixed" & tolerance.method != "auto") {
			stop("invalid argument for tolerance method - should be auto or fixed")
		}
	}

	if (verbose) {
		cat("ntrees resid. dev.","\n")
	}
	
	for (i in 1:n.folds) {

		model.mask <- selector != i  #used to fit model on majority of data
		pred.mask <- selector == i   #used to identify the with-held subset

		y.subset <- y.data[model.mask]
		x.subset <- x.data[model.mask, ,drop=FALSE]
		weight.subset <- site.weights[model.mask]

		if (!is.null(offset)) {
			offset.subset <- offset[model.mask] 
		} else {
			offset.subset <- NULL
		}

		model.list[[i]] <- eval(parse(text = gbm.call))

		fitted.values <- model.list[[i]]$fit  #predict.gbm(model.list[[i]], x.subset, type = "response", n.trees = n.trees)
		if (!is.null(offset)) {
			fitted.values <- fitted.values + offset[model.mask]
		}
		if (family == "bernoulli") {
			fitted.values <- exp(fitted.values)/(1 + exp(fitted.values))
		} else if (family == "poisson") {
			fitted.values <- exp(fitted.values)
		}
		
		pred.values[pred.mask] <- gbm::predict.gbm(model.list[[i]], x.data[pred.mask, ,drop=FALSE], n.trees = n.trees)
		
		if (!is.null(offset)) {
			pred.values[pred.mask] <- pred.values[pred.mask] + offset[pred.mask]
		}
		if (family == "bernoulli") {
			pred.values[pred.mask] <- exp(pred.values[pred.mask])/(1 + exp(pred.values[pred.mask]))
		} else if (family == "poisson") {
			pred.values[pred.mask] <- exp(pred.values[pred.mask])
		}

	# calc training deviance

		y_i <- y.subset
		u_i <- fitted.values
		weight.fitted <- site.weights[model.mask]
		training.loss.matrix[i,1] <- machisplin.calc.deviance(y_i, u_i, weight.fitted, family = family)

# calc holdout deviance

		y_i <- y.data[pred.mask]
		u_i <- pred.values[pred.mask]
		weight.preds <- site.weights[pred.mask]
		cv.loss.matrix[i,1] <- machisplin.calc.deviance(y_i, u_i, weight.preds, family = family)

	} # end of first loop

# now process until the change in mean deviance is =< tolerance or max.trees is exceeded

	delta.deviance <- 1

	cv.loss.values <- apply(cv.loss.matrix,2,mean)
	if (verbose) {
		cat(n.fitted,"  ",round(cv.loss.values,4),"\n")
	}
	if (!silent) {
		cat("now adding trees...","\n")
	}
	
	j <- 1

	while (delta.deviance > tolerance.test & n.fitted < max.trees) { 
		# beginning of inner loop

		# add a new column to the results matrice..

		training.loss.matrix <- cbind(training.loss.matrix,rep(0,n.folds))
		cv.loss.matrix <- cbind(cv.loss.matrix,rep(0,n.folds))

		n.fitted <- n.fitted + step.size
		trees.fitted <- c(trees.fitted,n.fitted)
  
		j <- j + 1
 
		for (i in 1:n.folds) {

			model.mask <- selector != i  #used to fit model on majority of data
			pred.mask <- selector == i   #used to identify the with-held subset

			y.subset <- y.data[model.mask]
			x.subset <- x.data[model.mask, ,drop=FALSE]
			weight.subset <- site.weights[model.mask]
			if (!is.null(offset)) {
				offset.subset <- offset[model.mask] 
			}
			model.list[[i]] <- gbm::gbm.more(model.list[[i]], weights = weight.subset, step.size)

			fitted.values <- model.list[[i]]$fit # predict.gbm(model.list[[i]],x.subset, type = "response", n.trees = n.fitted) 
			if (!is.null(offset)) {
				fitted.values <- fitted.values + offset[model.mask]
			}
			if (family == "bernoulli") {
				fitted.values <- exp(fitted.values)/(1 + exp(fitted.values))
			} else if (family == "poisson") {
				fitted.values <- exp(fitted.values)
			}
			pred.values[pred.mask] <- gbm::predict.gbm(model.list[[i]], x.data[pred.mask, ,drop=FALSE], n.trees = n.fitted)
			
			if (!is.null(offset)) {
				pred.values[pred.mask] <- pred.values[pred.mask] + offset[pred.mask]
			}
			
			if (family == "bernoulli") {
				pred.values[pred.mask] <- exp(pred.values[pred.mask])/(1 + exp(pred.values[pred.mask]))
			} else if (family == "poisson") {
				pred.values[pred.mask] <- exp(pred.values[pred.mask])
			}
# calculate training deviance

			y_i <- y.subset
			u_i <- fitted.values
			weight.fitted <- site.weights[model.mask]
			training.loss.matrix[i,j] <- machisplin.calc.deviance(y_i, u_i, weight.fitted, family = family)

# calc holdout deviance

			u_i <- pred.values[pred.mask]
			y_i <- y.data[pred.mask]
			weight.preds <- site.weights[pred.mask]
			cv.loss.matrix[i,j] <- machisplin.calc.deviance(y_i, u_i, weight.preds, family = family)

		}  # end of inner loop

		cv.loss.values <- apply(cv.loss.matrix,2,mean)

		if (j < 5) {
			if (cv.loss.values[j] > cv.loss.values[j-1]) {
				if (!silent) {
					cat("restart model with a smaller learning rate or smaller step size...")
				}
				return()
			}
		}

		if (j >= 20) {   #calculate stopping rule value
			test1 <- mean(cv.loss.values[(j-9):j])
			test2 <- mean(cv.loss.values[(j-19):(j-9)])
			delta.deviance <- test2 - test1
		}

		if (verbose) {
			cat(n.fitted," ",round(cv.loss.values[j],4),"\n") 
			flush.console()
		}
	} # end of while loop

# now begin process of calculating optimal number of trees

	training.loss.values <- apply(training.loss.matrix,2,mean)

	cv.loss.ses <- rep(0,length(cv.loss.values))
	cv.loss.ses <- sqrt(apply(cv.loss.matrix,2,var)) / sqrt(n.folds)

# find the target holdout deviance

	y.bar <- min(cv.loss.values) 


# identify the optimal number of trees 

	target.trees <- trees.fitted[match(TRUE, cv.loss.values == y.bar)]

# plot out the resulting curve of holdout deviance 
	if (plot.main) {

		y.min <- min(cv.loss.values - cv.loss.ses)  #je added multiplier 10/8/05
		y.max <- max(cv.loss.values + cv.loss.ses)  #je added multiplier 10/8/05 }

		if (plot.folds) {
			y.min <- min(cv.loss.matrix)
			y.max <- max(cv.loss.matrix) 
		}

		plot(trees.fitted, cv.loss.values, type = 'l', ylab = "holdout deviance", xlab = "no. of trees", ylim = c(y.min,y.max), ...)
		abline(h = y.bar, col = 2)

		lines(trees.fitted, cv.loss.values + cv.loss.ses, lty=2)  
		lines(trees.fitted, cv.loss.values - cv.loss.ses, lty=2)  

		if (plot.folds) {
			for (i in 1:n.folds) {
				lines(trees.fitted, cv.loss.matrix[i,],lty = 3)
			}
		}
		abline(v = target.trees, col=3)
		title(paste(sp.name,", d - ",tree.complexity,", lr - ",learning.rate, sep=""))
	}

# estimate the cv deviance and test statistics
# includes estimates of the standard error of the fitted values added 2nd may 2005

	cv.deviance.stats <- rep(0, n.folds)
	cv.roc.stats <- rep(0, n.folds)
	cv.cor.stats <- rep(0, n.folds)
	cv.calibration.stats <- matrix(0, ncol=5, nrow = n.folds)
	if (family == "bernoulli") {
		threshold.stats <- rep(0, n.folds)
	}
	fitted.matrix <- matrix(NA, nrow = n.cases, ncol = n.folds)  # used to calculate se's
	fold.fit <- rep(0, n.cases)

	for (i in 1:n.folds) {

		pred.mask <- selector == i   #used to identify the with-held subset
		model.mask <- selector != i  #used to fit model on majority of data

		fits <- gbm::predict.gbm(model.list[[i]], x.data[model.mask, ,drop=FALSE], n.trees = target.trees)
		if (!is.null(offset)) {
			fits <- fits + offset[model.mask]
		}
		if (family == "bernoulli") {
			fits <- exp(fits)/(1 + exp(fits))
		} else if (family == "poisson") {
			fits <- exp(fits)
		}
		fitted.matrix[model.mask,i] <- fits

		fits <- gbm::predict.gbm(model.list[[i]], x.data[pred.mask, ,drop=FALSE], n.trees = target.trees)
		if (!is.null(offset)) fits <- fits + offset[pred.mask]
		fold.fit[pred.mask] <- fits  # store the linear predictor values
		if (family == "bernoulli") {
			fits <- exp(fits)/(1 + exp(fits))
		} else if (family == "poisson") {
			fits <- exp(fits)
		}
		fitted.matrix[pred.mask,i] <- fits

		y_i <- y.data[pred.mask] 
		u_i <- fitted.matrix[pred.mask,i]  #pred.values[pred.mask]
		weight.preds <- site.weights[pred.mask]

		cv.deviance.stats[i] <- machisplin.calc.deviance(y_i, u_i, weight.preds, family = family)

		cv.cor.stats[i] <- cor(y_i,u_i)

		if (family == "bernoulli") {
			cv.roc.stats[i] <- .roc(y_i,u_i)
			cv.calibration.stats[i,] <- .calibration(y_i,u_i,"binomial")
			threshold.stats[i] <- approx(ppoints(u_i), sort(u_i,decreasing = T), prevalence)$y
		}

		if (family == "poisson") {
			cv.calibration.stats[i,] <- .calibration(y_i,u_i,"poisson")
		}
	}

	fitted.vars <- apply(fitted.matrix,1, var, na.rm = TRUE)

# now calculate the mean and se's for the folds

	cv.dev <- mean(cv.deviance.stats, na.rm = TRUE)
	cv.dev.se <- sqrt(var(cv.deviance.stats)) / sqrt(n.folds)

	cv.cor <- mean(cv.cor.stats, na.rm = TRUE)
	cv.cor.se <- sqrt(var(cv.cor.stats, use = "complete.obs")) / sqrt(n.folds)

	cv.roc <- 0.0
	cv.roc.se <- 0.0 

	if (family == "bernoulli") {
		cv.roc <- mean(cv.roc.stats,na.rm=TRUE)
		cv.roc.se <- sqrt(var(cv.roc.stats, use = "complete.obs")) / sqrt(n.folds)  
		cv.threshold <- mean(threshold.stats, na.rm = T)
		cv.threshold.se <- sqrt(var(threshold.stats, use = "complete.obs")) / sqrt(n.folds)  
	}
 
	cv.calibration <- 0.0
	cv.calibration.se <- 0.0

	if (family == "poisson" | family == "bernoulli") {
		cv.calibration <- apply(cv.calibration.stats,2,mean)
		cv.calibration.se <- apply(cv.calibration.stats,2,var)
		cv.calibration.se <- sqrt(cv.calibration.se) / sqrt(n.folds) 
	}

# fit the final model

	if (is.null(offset)) {
		gbm.call <- paste("gbm::gbm(y.data ~ .,data=x.data, n.trees = target.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = site.weights, distribution = as.character(family), var.monotone = var.monotone, verbose = FALSE)", sep="")
	} else {
		gbm.call <- paste("gbm::gbm(y.data ~ . + offset(offset),data=x.data, n.trees = target.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = site.weights, distribution = as.character(family), var.monotone = var.monotone,  verbose = FALSE)", sep="")
    } 

	if (!silent) {
		message("fitting final gbm model with a fixed number of ", target.trees, " trees for ", sp.name) 
	}
	gbm.object <- eval(parse(text = gbm.call))

	best.trees <- target.trees

#extract fitted values and summary table
  
	gbm.summary <- summary(gbm.object,n.trees = target.trees, plotit = FALSE)

	fits <- gbm::predict.gbm(gbm.object, x.data, n.trees = target.trees)
	if (!is.null(offset)) fits <- fits + offset
	if (family == "bernoulli") {
		fits <- exp(fits)/(1 + exp(fits))
	} else if (family == "poisson") {
		fits <- exp(fits)
	}
	fitted.values <- fits

	y_i <- y.data
	u_i <- fitted.values
	resid.deviance <- machisplin.calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)

	self.cor <- cor(y_i,u_i)
	self.calibration <- 0.0
	self.roc <- 0.0

	if (family == "bernoulli") {  # do this manually as we need the residuals
		deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
		residuals <- sqrt(abs(deviance.contribs * 2))
		residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
		self.roc <- .roc(y_i,u_i)
		self.calibration <- .calibration(y_i,u_i,"binomial")
	}

	if (family == "poisson") {   # do this manually as we need the residuals
		deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
		residuals <- sqrt(abs(deviance.contribs * 2))
		residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
		self.calibration <- .calibration(y_i,u_i,"poisson")
	}

	if (family == "gaussian" | family == "laplace") {
		residuals <- y_i - u_i
	}

	mean.resid.deviance <- resid.deviance/n.cases

	z2 <- Sys.time()
	elapsed.time.minutes <- round(as.numeric(z2 - z1)/ 60, 2)  #calculate the total elapsed time

	if (verbose) {
		cat("\n")
		cat("mean total deviance =", round(mean.total.deviance,3),"\n")
		cat("mean residual deviance =", round(mean.resid.deviance,3),"\n","\n")
		cat("estimated cv deviance =", round(cv.dev,3),"; se =", round(cv.dev.se,3),"\n","\n")
		cat("training data correlation =",round(self.cor,3),"\n")
		cat("cv correlation = ",round(cv.cor,3),"; se =",round(cv.cor.se,3),"\n","\n")
		if (family == "bernoulli") {
			cat("training data AUC score =",round(self.roc,3),"\n")
			cat("cv AUC score =",round(cv.roc,3),"; se =",round(cv.roc.se,3),"\n","\n")
		}
		cat("elapsed time - ",round(elapsed.time.minutes,2),"minutes","\n")
	}

	if (n.fitted == max.trees & !silent) {
		cat("\n","########### warning ##########","\n","\n")
		cat("maximum tree limit reached - results may not be optimal","\n")
		cat("  - refit with faster learning rate or increase maximum number of trees","\n") 
	}
    
# now assemble data to be returned

	gbm.detail <- list(dataframe = data, gbm.x = gbm.x, predictor.names = names(x.data), 
    gbm.y = gbm.y, response.name = sp.name, offset = offset, family = family, tree.complexity = tree.complexity, 
    learning.rate = learning.rate, bag.fraction = bag.fraction, cv.folds = n.folds, 
    prev.stratification = prev.stratify, max.fitted = n.fitted, n.trees = target.trees, 
    best.trees = target.trees, train.fraction = 1.0, tolerance.method = tolerance.method, 
    tolerance = tolerance, var.monotone = var.monotone, date = date(), 
    elapsed.time.minutes = elapsed.time.minutes)

	training.stats <- list(null = total.deviance, mean.null = mean.total.deviance, 
    resid = resid.deviance, mean.resid = mean.resid.deviance, correlation = self.cor, 
    discrimination = self.roc, calibration = self.calibration)

	cv.stats <- list(deviance.mean = cv.dev, deviance.se = cv.dev.se, 
    correlation.mean = cv.cor, correlation.se = cv.cor.se,
    discrimination.mean = cv.roc, discrimination.se = cv.roc.se,
    calibration.mean = cv.calibration, calibration.se = cv.calibration.se)

	if (family == "bernoulli") {
		cv.stats$cv.threshold <- cv.threshold
		cv.stats$cv.threshold.se <- cv.threshold.se
	}

#  rm(x.data,y.data, envir = globalenv())           #finally, clean up the temporary dataframes

# and assemble results for return

	gbm.object$gbm.call <- gbm.detail
	gbm.object$fitted <- fitted.values
	gbm.object$fitted.vars <- fitted.vars
	gbm.object$residuals <- residuals
	gbm.object$contributions <- gbm.summary
	gbm.object$self.statistics <- training.stats
	gbm.object$cv.statistics <- cv.stats
	gbm.object$weights <- site.weights
	gbm.object$trees.fitted <- trees.fitted
	gbm.object$training.loss.values <- training.loss.values
	gbm.object$cv.values <- cv.loss.values
	gbm.object$cv.loss.ses <- cv.loss.ses
	gbm.object$cv.loss.matrix <- cv.loss.matrix
	gbm.object$cv.roc.matrix <- cv.roc.stats

	if (keep.fold.models) {
		gbm.object$fold.models <- model.list
	} else {
		gbm.object$fold.models <- NULL
	}

	if (keep.fold.vector) {
		gbm.object$fold.vector <- selector
	} else {
		gbm.object$fold.vector <- NULL
	}
	if (keep.fold.fit) {
		gbm.object$fold.fit <- fold.fit
	} else {
		gbm.object$fold.fit <- NULL
	}
	
	return(gbm.object)  
}

# j. leathwick/j. elith
#
# version 2.1 - 5th Sept 2005
#
# function to calculate deviance given two vectors of raw and fitted values
# requires a family argument which is set to binomial by default
#
#

machisplin.calc.deviance <-  function(obs, pred, weights = rep(1,length(obs)), family="binomial", calc.mean = TRUE) {

if (length(obs) != length(pred)) {   stop("observations and predictions must be of equal length") }

y_i <- obs
u_i <- pred
 
family = tolower(family)
 
if (family == "binomial" | family == "bernoulli") {
 
   deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
   deviance <- -2 * sum(deviance.contribs * weights)

} else if (family == "poisson") {

    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    deviance <- 2 * sum(deviance.contribs * weights)

} else if (family == "laplace") {

    deviance <- sum(abs(y_i - u_i))
	
} else if (family == "gaussian") {

    deviance <- sum((y_i - u_i) * (y_i - u_i))
	
} else {
	stop('unknown family, should be one of: "binomial", "bernoulli", "poisson", "laplace", "gaussian"')
}

if (calc.mean) deviance <- deviance/length(obs)

return(deviance)

}