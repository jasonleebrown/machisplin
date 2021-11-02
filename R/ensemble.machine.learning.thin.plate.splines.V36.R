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
#' @import raster
#' @import dismo
#' @import snow
#' @import snowfall
#' @import maptools
#' @import rgdal
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
#' library(raster)
#' 
#' # Import a csv as shapefile:
#' Mydata<-read.delim("sampling.csv", sep=",", h=T)
#' 
#' #load rasters to use as high resolution co-variates for downscaling
#' ALT = raster("SRTM30m.tif")
#' SLOPE = raster("ln_slope.tif")
#' ASPECT = raster("aspect.tif")
#' GEOMORPH = raster("geomorphons.tif")
#' TWI = raster("TWI.tif")
#' 
#' # function input: raster brick of covarites
#' raster_stack<-stack(ALT,SLOPE,TWI,GEOMORPH, ASPECT)
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

machisplin.mltps<-function(int.values, covar.ras, n.cores=1, tps=TRUE, smooth.outputs.only=FALSE){
f <- list()

i.lyrs<-ncol(int.values)
# n of iterpolations - subtrack x and y
n.spln<-i.lyrs-2

# n of iterpolations-add two layers for lat and long
n.covars<-nlayers(covar.ras)+2

## create a raster for lat and long and any other covariate, add to raster stack
## create rasters of LAT and LONG for downscalling, using ALT as template for dimensions
r<-covar.ras[[1]]
LAT <- LONG <- r
xy <- coordinates(r)
LONG[] <- xy[, 1]
LAT[] <- xy[, 2]

#########################################
#load raster for statistical downscaling (a raster brick saved in the grid format and countaining all topographic data)
rast.names<-names(covar.ras)
rast.names.all<-c(rast.names,"LONG","LAT")
rast_stack<-stack(covar.ras,LONG, LAT)
n.names<-(1:nlayers(rast_stack))
for(i in n.names){
names(rast_stack)[[i]]<- rast.names.all[i]}

#extract values to points from rasters
RAST_VAL<-data.frame(extract(rast_stack, int.values[1:2]))
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
MyShp<-SpatialPointsDataFrame(Mydata[,1:2],Mydata)

#specify coordinate system
crs(MyShp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
 

#######################
#make list to loop over 
#the column 3 to 4 in temp_dat are the temperature variable 
#and the column 5 to xx are the topographic variables 
dat_tps<-list()
n.bios<-c(1:n.spln)# number of climate variables, here two
for (i in n.bios){
    dat_tps[[i]]<-MyShp[,c((i+2), c((i.lyrs+1):(n.covars+i.lyrs)))] #2 is b/c 3 and 4 col are temp data thus i+2=3 & 4 columns
  }
#str(dat_tps)

#rename each dataset within the list
name.dat<-c("resp",rast.names.all)
dat_tps<- lapply(seq(dat_tps), function(i) {
        y <- data.frame(dat_tps[[i]])
        names(y) <- c(name.dat)
        return(y)
    })
#str(dat_tps)

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
			l <- list()
			nfolds <- 10
			kfolds <- kfold(dat_tps[[i]], nfolds)
						
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
			mod.brt.tps.elev<- gbm.step(data=train, gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 25, learning.rate = 0.01, bag.fraction = 0.5, plot.main = FALSE)
		    mod.rf.tps.elev<- randomForest(mod.form, data = train)
            mod.nn.tps.elev<-nnet(mod.form, data = trainNN, size=10, linout=TRUE, maxit=10000)
			mod.mars.tps.elev<-earth(mod.form, data = train, nfold=10)
            mod.svm.tps.elev<- ksvm(mod.form, data = train)
			mod.gam.tps.elev<- gam(mod.form, data = train)
			
			
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
				
				OptX<-optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4], x[5], x[6],mfit.brt.full,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.rf.full,mfit.svm.full), lower=0, upper=1)
				
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
				
				OptX<-optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4],mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.svm.full), lower=0, upper=1)
				
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
			    mod.nn.tps.FINAL<-nnet(mod.form, data = trainNN.f, size=10, linout=TRUE, maxit=10000)
			    #store variable importance
                l$var.imp$nn<-garson(mod.nn.tps.FINAL, bar_plot=F)
				
			    #create raster pred
			    if(is.null(pred.elev)==FALSE){pred.elev.nn<-predict(rast_stack, mod.nn.tps.FINAL)
			    pred.nn<-pred.elev.nn*max2.resp.f
		        pred.elev.2<-pred.nn+min.resp.f
			    pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
			    if(is.null(pred.elev)==TRUE){pred.elev.nn<-predict(rast_stack, mod.nn.tps.elev)
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
				mod.brt.tps.FINAL<- gbm.step(data=dat_tps[[i]], gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5, plot.main = FALSE)
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
				mod.rf.tps.FINAL<- randomForest(mod.form, data = dat_tps[[i]],importance = TRUE)
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
				mod.MARS.tps.FINAL<-earth(mod.form,  data = dat_tps[[i]], nfold=10)
				#store variable importance
				l$var.imp$mars<-evimp(mod.MARS.tps.FINAL)
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
				mod.SVM.tps.FINAL<-ksvm(mod.form, data=dat_tps[[i]])
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
				mod.GAM.tps.FINAL<-gam(mod.form, data=dat_tps[[i]])
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
						c <- as(extent(new.df2[[h]][1], new.df2[[h]][2], new.df2[[h]][3], new.df2[[h]][4]), 'SpatialPolygons')
						crs(b) <- crs(rast_stack)
						crs(c) <- crs(rast_stack)
						rb <- crop(rast_stack, b)
						#sub sample residuals
						RAST_TPS<-data.frame(extract(rb[[1]], Full.cords))
						#str(RAST_VAL)
						#merge sampled data to input
						MyTPSdata<-cbind(res.FINAL,Full.cords,RAST_TPS)
						#remove NAs
						MyTPSdata<-MyTPSdata[complete.cases(MyTPSdata),][1:3]
						#fit thin plate spline of residuals
						mod.tps.elev<-Tps(MyTPSdata[2:3], MyTPSdata[1])#columns of lat and long
						#use TPS to interpolate residuals
						TPS_name<-paste0("TILE_",h)
						#use TPS to interpolate residuals
						pred_TPS_elev<-interpolate(rb, mod.tps.elev)
						names(pred_TPS_elev)<- TPS_name
						pred_TPS_elev <- crop(pred_TPS_elev, c)
						pred_TPS_elev <- extend(pred_TPS_elev,rast_stack)
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
				l$residuals<- cbind((dat_tps[[i]]$resp - f.actual),dat_tps[[i]]$LONG,dat_tps[[i]]$LAT)
				colnames(l$residuals)<-c("residuals","long","lat")
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
				if(i==1){l$n.layers<-length(out.names)}
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
	if(n.spln>1){i<-seq(1,n.spln)} else {i<-1}# length = number of climate variables
			##################################################################################################
			############################# part 1 evaluate best ensemble of models ##############################
			##################################################################################################
            #perform K-fold cross val
			l <- list()
			nfolds <- 10
			kfolds <- kfold(dat_tps[[i]], nfolds)
						
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
			mod.brt.tps.elev<- gbm.step(data=train, gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 25, learning.rate = 0.01, bag.fraction = 0.5, plot.main = FALSE)
		    mod.rf.tps.elev<- randomForest(mod.form, data = train)
            mod.nn.tps.elev<-nnet(mod.form, data = trainNN, size=10, linout=TRUE, maxit=10000)
			mod.mars.tps.elev<-earth(mod.form, data = train, nfold=10)
            mod.svm.tps.elev<- ksvm(mod.form, data = train)
			mod.gam.tps.elev<- gam(mod.form, data = train)
			
			
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
				
				OptX<-optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4], x[5], x[6],mfit.brt.full,mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.rf.full,mfit.svm.full), lower=0, upper=1)
				
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
				
				OptX<-optimx(par, function(x) machisplin.optimx.internal(x[1], x[2], x[3], x[4],mfit.gam.full,mfit.nn.full,mfit.mars.full,mfit.svm.full), lower=0, upper=1)
				
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
			    mod.nn.tps.FINAL<-nnet(mod.form, data = trainNN.f, size=10, linout=TRUE, maxit=10000)
			    #store variable importance
                l$var.imp$nn<-garson(mod.nn.tps.FINAL, bar_plot=F)
				
			    #create raster pred
			    if(is.null(pred.elev)==FALSE){pred.elev.nn<-predict(rast_stack, mod.nn.tps.FINAL)
			    pred.nn<-pred.elev.nn*max2.resp.f
		        pred.elev.2<-pred.nn+min.resp.f
			    pred.elev<-(pred.elev+(pred.elev.2*OptX.mfit.wt[[iter.mod]]))}
			    if(is.null(pred.elev)==TRUE){pred.elev.nn<-predict(rast_stack, mod.nn.tps.elev)
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
				mod.brt.tps.FINAL<- gbm.step(data=dat_tps[[i]], gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5, plot.main = FALSE)
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
				mod.rf.tps.FINAL<- randomForest(mod.form, data = dat_tps[[i]],importance = TRUE)
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
				mod.MARS.tps.FINAL<-earth(mod.form,  data = dat_tps[[i]], nfold=10)
				#store variable importance
				l$var.imp$mars<-evimp(mod.MARS.tps.FINAL)
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
				mod.SVM.tps.FINAL<-ksvm(mod.form, data=dat_tps[[i]])
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
				mod.GAM.tps.FINAL<-gam(mod.form, data=dat_tps[[i]])
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
						c <- as(extent(new.df2[[h]][1], new.df2[[h]][2], new.df2[[h]][3], new.df2[[h]][4]), 'SpatialPolygons')
						crs(b) <- crs(rast_stack)
						crs(c) <- crs(rast_stack)
						rb <- crop(rast_stack, b)
						#sub sample residuals
						RAST_TPS<-data.frame(extract(rb[[1]], Full.cords))
						#str(RAST_VAL)
						#merge sampled data to input
						MyTPSdata<-cbind(res.FINAL,Full.cords,RAST_TPS)
						#remove NAs
						MyTPSdata<-MyTPSdata[complete.cases(MyTPSdata),][1:3]
						#fit thin plate spline of residuals
						mod.tps.elev<-Tps(MyTPSdata[2:3], MyTPSdata[1])#columns of lat and long
						#use TPS to interpolate residuals
						TPS_name<-paste0("TILE_",h)
						#use TPS to interpolate residuals
						pred_TPS_elev<-interpolate(rb, mod.tps.elev)
						names(pred_TPS_elev)<- TPS_name
						pred_TPS_elev <- crop(pred_TPS_elev, c)
						pred_TPS_elev <- extend(pred_TPS_elev,rast_stack)
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
				l$residuals<- cbind((dat_tps[[i]]$resp - f.actual),dat_tps[[i]]$LONG,dat_tps[[i]]$LAT)
				colnames(l$residuals)<-c("residuals","long","lat")
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
				if(i==1){l$n.layers<-length(out.names)}
				l$summary<-sum.val
				#save TPS if R2 is better, else save only model
				names(pred.elev)<- out.names[i]
				l$final<-pred.elev
				}
		l$var.imp
		return(l)
	}
}

##################################################################################################
##################################################################################################
#' Write MACHISPLIN geotiff
#' @param mttps.in an output from the 'machisplin.mltps' function  
#' @param out.names a vector corresponding to the output raster names.   If 'null' it will write intial input names.  
#' @return This function outputs performance values, the algorithm(s) used, and rasters for use in GIS output from the 'machisplin.mltps' function
#' @export
#' @import raster
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
		if(i ==1){rast.O<-mltps.in[[i]]$final} else {rast.O<-stack(rast.O, mltps.in[[i]]$final)}
		}
	#write rasters
	if(is.null(out.names)==TRUE){writeRaster(stack(rast.O), names(rast.O), bylayer=TRUE, format='GTiff',overwrite=overwrite)}
	if(is.null(out.names)==FALSE){writeRaster(stack(rast.O), names(out.names), bylayer=TRUE, format='GTiff',overwrite=overwrite)}

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
#' @import raster
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
		print(interp.rast[[i]]$var.imp)
		sink()		
		}
}
 # returns output to the console
##################################################################################################
##################################################################################################
#' Write MACHISPLIN model residuals
#' @param mttps.in an output from the 'machisplin.mltps' function  
#' @param out.names a vector corresponding to the output raster names.   If 'null' it will write intial input names.  
#' @return This function outputs performance values, the algorithm(s) used, and rasters for use in GIS output from the 'machisplin.mltps' function
#' @export
#' @import raster
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
#import: nnet and optimx
"
library(raster)
library(gstat)
library(MASS)
library(nlme)
library(fields)
library(dismo)
library(randomForest)
library(earth)
library(kernlab)
library(mgcv)
library(snowfall)
library(snow)
library(maptools)
library(rgdal)
library(spdep)
library(nnet)
library(optimx)
"
