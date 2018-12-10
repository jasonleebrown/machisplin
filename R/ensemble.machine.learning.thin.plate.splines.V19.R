#add features to deal with empty grid for tilted TPS
##################################################################################################
##################################################################################################
#' Machine Learning Ensemble & Thin-Plate-Spline Interpolation
#' @param imported table for analysis  
#' @return This tool remove NAs and converts all input data to numeric values for analysis in Humboldt. 
#' @export
#' @examples
#' library(machisplin)
#' 
#' library(raster)
#' library(gstat)
#' library(maptools)
#' library(sp)
#' library(rgdal)
#' library(spdep)
#' library(nlme)
#' library(MASS)
#' library(fields)
#' library(snowfall)
#' library(dismo)
#' library(randomForest)
#' library(earth)
#' library(kernlab)
#' library(mgcv)
#' 
#' # Import a csv as shapefile:
#' Mydata<-read.delim("sampling.csv", sep=",", h=T)
#' 
#' #load rasters
#' ALT = raster("SRTM30m.tif")
#' CURV_V = raster("cs_curv.tif")
#' CURV_H = raster("long_cur.tif")
#' SLOPE = raster("ln_slope.tif")
#' REL_ALT= raster("relative_elevation500m.tif")
#' TWI = raster("TWI.tif")
#' 
#' # function input: raster brick of covarites
#' raster_stack<-stack(ALT,REL_ALT,SLOPE,TWI,CURV_V,CURV_H)

#' 
#' #run an ensemble machine learning thin plate spline 
#' interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2)
#' 


machisplin.mltps<-function(int.values, covar.ras, n.cores=1){
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
RAST_VAL<-data.frame(extract(rast_stack, Mydata[1:2]))
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
out.names<-names(Mydata[,nm.spln.out])

#extract names for model formula
xnam <- name.dat[2:(n.covars+1)]
mod.form<- as.formula(paste("resp ~ ", paste(xnam, collapse= "+")))

########################################################################
#########################################################################
#function to pass to snowfall for parallel computing 
if(n.spln>1){i<-seq(1,n.spln)} else {i<-1}# length = number of climate variables
  myLapply_elevOut<-function(i){
			##################################################################################################
			############################# part 1 evaluate best ensemble of models ##############################
			##################################################################################################
            #perform K-fold cross val
			l <- list()
			nfolds <- 10
			kfolds <- kfold(dat_tps[[i]], nfolds)
			
			mfit.brt <- mfit.rf <- mfit.l <- mfit.mars<- mfit.svm <- mfit.gam <- mfit.br<-mfit.bl<-mfit.bm<-mfit.lr<-mfit.lm<-mfit.rm<- mfit.bg<-mfit.rg<-mfit.lg<-mfit.mg<-mfit.sg<-mfit.bv<-mfit.rv<-mfit.lv<-mfit.mv<-mfit.rmv<-mfit.rbv<-mfit.rlv<-mfit.rmg<-mfit.rbg<-mfit.rlg<-mfit.rgv<-mfit.blv<-mfit.bmv<-mfit.blg<-mfit.bmg<-mfit.mlg<-mfit.mlv<-mfit.mvg<-mfit.lvg<-mfit.bvg<-mfit.rml<-mfit.rmb<-mfit.rbl<-mfit.mbl<-mfit.rmbg<-mfit.rmbv<-mfit.rmbl<-mfit.vmbg<-mfit.lmbv<-mfit.gvlb<-mfit.gmlb<-mfit.grlb<-mfit.gvrb<-mfit.gvrm<-mfit.mvrl<-mfit.gmlv<-mfit.grlm<-mfit.grlv<-mfit.rlbv<-mfit.glmrv<-mfit.blmrv<-mfit.bgmrv<-mfit.bglrv<-mfit.bglmv<-mfit.bglmr<-mfit.bglmrv<- rep(NA, nfolds)
			
			for (v in 1:nfolds) {
			#train with 90% of data is <4000 or 10% if >4000 total rows
			if (nrow(dat_tps[[i]])>4000){
				train <- dat_tps[[i]][kfolds==v,]
				test <- dat_tps[[i]][kfolds!=v,]
			} else {train <- dat_tps[[i]][kfolds!=v,]
				test <- dat_tps[[i]][kfolds==v,]}
			####train models
            mod.brt.tps.elev<- gbm.step(data=train, gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 25, learning.rate = 0.01, bag.fraction = 0.5, plot.main = FALSE)
		    mod.rf.tps.elev<- randomForest(mod.form, data = train)
            mod.lm.tps.elev<-stepAIC(lm(mod.form, data = train))
			mod.mars.tps.elev<-earth(mod.form, data = train, nfold=10, ncross=30, varmod.method="lm")
            mod.svm.tps.elev<- ksvm(mod.form, data = train)
			mod.gam.tps.elev<- gam(mod.form, data = train)
			
			#extract residuals and calculate residual sum of squares
			#BRT
			pred.brt.obs <- predict(mod.brt.tps.elev, test, n.trees=mod.brt.tps.elev$gbm.call$best.trees, type="response")
            res.brt.elev<-test[,1]-pred.brt.obs
			mfit.brt[[v]]<-sum(res.brt.elev ^ 2)
			#RF
            pred.rf.obs<-predict(mod.rf.tps.elev, test, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
			res.rf.elev<-test[,1]-pred.rf.obs
			mfit.rf[[v]]<-sum(res.rf.elev ^ 2)
			#LM
			res.reg.elev<-test[,1]-predict(mod.lm.tps.elev, test)
			mfit.l[[v]]<-sum(res.reg.elev ^ 2)
			#MAR
			pred.mars.test<-predict(mod.mars.tps.elev,test)
			res.mars.elev<-as.vector(test[,1]-pred.mars.test) 
			mfit.mars[[v]]<-sum(res.mars.elev ^ 2)
			#SVM
			pred.svm.test<-predict(mod.svm.tps.elev,test)
			res.svm.elev<-test[,1]-pred.svm.test 
			mfit.svm[[v]]<-sum(res.svm.elev ^ 2)
			#GAM
			pred.gam.test<-predict(mod.gam.tps.elev,test)
			res.gam.elev<-as.vector(test[,1]-pred.gam.test)
			mfit.gam[[v]]<-sum(res.gam.elev ^ 2)
			
			#estimate ensemble model fits and residual sum of squares
	        mfit.br[[v]]<-sum(((res.brt.elev+res.rf.elev)/2)^ 2)
            mfit.bl[[v]]<-sum(((res.brt.elev+res.reg.elev)/2)^ 2)
			mfit.bm[[v]]<-sum(((res.brt.elev+res.mars.elev)/2)^ 2)
			mfit.lr[[v]]<-sum(((res.rf.elev+res.reg.elev)/2)^ 2)
			mfit.lm[[v]]<-sum(((res.mars.elev+res.reg.elev)/2)^ 2)
			mfit.rm[[v]]<-sum(((res.mars.elev+res.rf.elev)/2)^ 2)		
			mfit.bg[[v]]<-sum(((res.gam.elev+res.brt.elev)/2)^ 2)
			mfit.rg[[v]]<-sum(((res.gam.elev+res.rf.elev)/2)^ 2)
			mfit.lg[[v]]<-sum(((res.gam.elev+res.reg.elev)/2)^ 2)
			mfit.mg[[v]]<-sum(((res.gam.elev+res.mars.elev)/2)^ 2)
			mfit.sg[[v]]<-sum(((res.gam.elev+res.svm.elev)/2)^ 2)
			mfit.bv[[v]]<-sum(((res.svm.elev+res.brt.elev)/2)^ 2)
			mfit.rv[[v]]<-sum(((res.svm.elev+res.rf.elev)/2)^ 2)
			mfit.lv[[v]]<-sum(((res.svm.elev+res.reg.elev)/2)^ 2)
			mfit.mv[[v]]<-sum(((res.svm.elev+res.mars.elev)/2)^ 2)
			
			mfit.rmv[[v]]<-sum(((res.mars.elev+res.rf.elev+res.svm.elev)/3)^ 2)
			mfit.rbv[[v]]<-sum(((res.rf.elev+res.brt.elev+res.svm.elev)/3)^ 2)
			mfit.rlv[[v]]<-sum(((res.rf.elev+res.reg.elev+res.svm.elev)/3)^ 2)
			mfit.rgv[[v]]<-sum(((res.rf.elev+res.gam.elev+res.svm.elev)/3)^ 2)
			mfit.rmg[[v]]<-sum(((res.mars.elev+res.rf.elev+res.gam.elev)/3)^ 2)
			mfit.rbg[[v]]<-sum(((res.brt.elev+res.rf.elev+res.gam.elev)/3)^ 2)
			mfit.rlg[[v]]<-sum(((res.rf.elev+res.reg.elev+res.gam.elev)/3)^ 2)
			mfit.blv[[v]]<-sum(((res.brt.elev+res.reg.elev+res.svm.elev)/3)^ 2)
			mfit.bmv[[v]]<-sum(((res.brt.elev+res.mars.elev+res.svm.elev)/3)^ 2)
			mfit.blg[[v]]<-sum(((res.brt.elev+res.reg.elev+res.gam.elev)/3)^ 2)
			mfit.bmg[[v]]<-sum(((res.brt.elev+res.mars.elev+res.gam.elev)/3)^ 2)
			mfit.mlg[[v]]<-sum(((res.mars.elev+res.gam.elev+res.reg.elev)/3)^ 2)
			mfit.mlv[[v]]<-sum(((res.mars.elev+res.svm.elev+res.reg.elev)/3)^ 2)
			mfit.mvg[[v]]<-sum(((res.mars.elev+res.svm.elev+res.gam.elev)/3)^ 2)
			mfit.lvg[[v]]<-sum(((res.reg.elev+res.svm.elev+res.gam.elev)/3)^ 2)
  			mfit.bvg[[v]]<-sum(((res.brt.elev+res.svm.elev+res.gam.elev)/3)^ 2)
     		mfit.rml[[v]]<-sum(((res.mars.elev+res.rf.elev+res.reg.elev)/3)^ 2)
			mfit.rmb[[v]]<-sum(((res.mars.elev+res.rf.elev+res.brt.elev)/3)^ 2)
			mfit.rbl[[v]]<-sum(((res.brt.elev+res.rf.elev+res.reg.elev)/3)^ 2)
			mfit.mbl[[v]]<-sum(((res.mars.elev+res.reg.elev+res.brt.elev)/3)^ 2)
			mfit.rmbg[[v]]<-sum(((res.rf.elev+res.mars.elev+res.brt.elev+res.gam.elev)/4)^ 2)
			mfit.rmbv[[v]]<-sum(((res.mars.elev+res.rf.elev+res.brt.elev+res.svm.elev)/4)^ 2)
			mfit.rmbl[[v]]<-sum(((res.mars.elev+res.reg.elev+res.brt.elev+res.rf.elev)/4)^ 2)
			mfit.vmbg[[v]]<-sum(((res.mars.elev+res.gam.elev+res.brt.elev+res.svm.elev)/4)^ 2)
			mfit.lmbv[[v]]<-sum(((res.mars.elev+res.reg.elev+res.brt.elev+res.svm.elev)/4)^ 2)
			mfit.gvlb[[v]]<-sum(((res.gam.elev+res.reg.elev+res.brt.elev+res.svm.elev)/4)^ 2)
			mfit.gmlb[[v]]<-sum(((res.mars.elev+res.reg.elev+res.brt.elev+res.gam.elev)/4)^ 2)
        	mfit.grlb[[v]]<-sum(((res.gam.elev+res.reg.elev+res.brt.elev+res.rf.elev)/4)^ 2)
        	mfit.gvrb[[v]]<-sum(((res.gam.elev+res.svm.elev+res.brt.elev+res.rf.elev)/4)^ 2)
			mfit.gvrm[[v]]<-sum(((res.mars.elev+res.gam.elev+res.svm.elev+res.rf.elev)/4)^ 2)
			mfit.mvrl[[v]]<-sum(((res.mars.elev+res.reg.elev+res.svm.elev+res.rf.elev)/4)^ 2)
			mfit.gmlv[[v]]<-sum(((res.mars.elev+res.reg.elev+res.gam.elev+res.svm.elev)/4)^ 2)
			mfit.grlm[[v]]<-sum(((res.mars.elev+res.reg.elev+res.gam.elev+res.rf.elev)/4)^ 2)
			mfit.grlv[[v]]<-sum(((res.gam.elev+res.reg.elev+res.svm.elev+res.rf.elev)/4)^ 2)
			mfit.rlbv[[v]]<-sum(((res.svm.elev+res.reg.elev+res.brt.elev+res.rf.elev)/4)^ 2)
			
			mfit.glmrv[[v]]<-sum(((res.gam.elev+res.reg.elev+res.mars.elev+res.rf.elev+res.svm.elev)/5)^ 2)
			mfit.blmrv[[v]]<-sum(((res.brt.elev+res.reg.elev+res.mars.elev+res.rf.elev+res.svm.elev)/5)^ 2)
			mfit.bgmrv[[v]]<-sum(((res.brt.elev+res.gam.elev+res.mars.elev+res.rf.elev+res.svm.elev)/5)^ 2)
			mfit.bglrv[[v]]<-sum(((res.brt.elev+res.gam.elev+res.reg.elev+res.rf.elev+res.svm.elev)/5)^ 2)
			mfit.bglmv[[v]]<-sum(((res.brt.elev+res.gam.elev+res.reg.elev+res.mars.elev+res.svm.elev)/5)^ 2)
			mfit.bglmr[[v]]<-sum(((res.brt.elev+res.gam.elev+res.reg.elev+res.mars.elev+res.rf.elev)/5)^ 2)
		
			mfit.bglmrv[[v]]<-sum(((res.brt.elev+res.gam.elev+res.reg.elev+res.mars.elev+res.rf.elev+res.svm.elev)/6)^ 2)
			}
			
			#small function to calculate AICc of Ensembles, NEZ is num of models, datz= dat_tps, nfitz=mfit...
			machisplin.aicc<-function(mfitz,datz,nez=1){
			mfitA<-mean(mfitz)
			mfitB<-nrow(datz)*log(mfitA/nrow(datz))+2*nez
			mfitC<-mfitB+(2*nez*(nez+1))/(nrow(dat_tps[[i]])-nez-1)
			return(mfitC)
			}
			#average k-folds
			mfit.brt <-machisplin.aicc(mfit.brt,dat_tps[[i]],1)
			mfit.rf <- machisplin.aicc(mfit.rf,dat_tps[[i]],1)
			mfit.l <- machisplin.aicc(mfit.lm,dat_tps[[i]],1)
			mfit.mars<- machisplin.aicc(mfit.mars,dat_tps[[i]],1)
			mfit.svm <- machisplin.aicc(mfit.svm,dat_tps[[i]],1)
			mfit.gam <- machisplin.aicc(mfit.gam,dat_tps[[i]],1)
			mfit.br<-machisplin.aicc(mfit.br,dat_tps[[i]],2)
			mfit.bl<-machisplin.aicc(mfit.bl,dat_tps[[i]],2)
			mfit.bm<-machisplin.aicc(mfit.bm,dat_tps[[i]],2)
			mfit.lr<-machisplin.aicc(mfit.lr,dat_tps[[i]],2)
			mfit.lm<-machisplin.aicc(mfit.lm,dat_tps[[i]],2)
			mfit.rm<- machisplin.aicc(mfit.rm,dat_tps[[i]],2)
			mfit.bg<-machisplin.aicc(mfit.bg,dat_tps[[i]],2)
			mfit.rg<-machisplin.aicc(mfit.rg,dat_tps[[i]],2)
			mfit.lg<-machisplin.aicc(mfit.lg,dat_tps[[i]],2)
			mfit.mg<-machisplin.aicc(mfit.mg,dat_tps[[i]],2)
			mfit.sg<-machisplin.aicc(mfit.sg,dat_tps[[i]],2)
			mfit.bv<-machisplin.aicc(mfit.bv,dat_tps[[i]],2)
			mfit.rv<-machisplin.aicc(mfit.rv,dat_tps[[i]],2)
			mfit.lv<-machisplin.aicc(mfit.lv,dat_tps[[i]],2)
			mfit.mv<-machisplin.aicc(mfit.mv,dat_tps[[i]],2)
			mfit.rmv<-machisplin.aicc(mfit.rmv,dat_tps[[i]],3)
			mfit.rbv<-machisplin.aicc(mfit.rbv,dat_tps[[i]],3)
			mfit.rlv<-machisplin.aicc(mfit.rlv,dat_tps[[i]],3)
			mfit.rmg<-machisplin.aicc(mfit.rmg,dat_tps[[i]],3)
			mfit.rbg<-machisplin.aicc(mfit.rbg,dat_tps[[i]],3)
			mfit.rlg<-machisplin.aicc(mfit.rlg,dat_tps[[i]],3)
			mfit.rgv<-machisplin.aicc(mfit.rgv,dat_tps[[i]],3)
			mfit.blv<-machisplin.aicc(mfit.blv,dat_tps[[i]],3)
			mfit.bmv<-machisplin.aicc(mfit.bmv,dat_tps[[i]],3)
			mfit.blg<-machisplin.aicc(mfit.blg,dat_tps[[i]],3)
			mfit.bmg<-machisplin.aicc(mfit.bmg,dat_tps[[i]],3)
			mfit.mlg<-machisplin.aicc(mfit.mlg,dat_tps[[i]],3)
			mfit.mlv<-machisplin.aicc(mfit.mlv,dat_tps[[i]],3)
			mfit.mvg<-machisplin.aicc(mfit.mvg,dat_tps[[i]],3)
			mfit.lvg<-machisplin.aicc(mfit.lvg,dat_tps[[i]],3)
			mfit.bvg<-machisplin.aicc(mfit.bvg,dat_tps[[i]],3)
			mfit.rml<-machisplin.aicc(mfit.rml,dat_tps[[i]],3)
			mfit.rmb<-machisplin.aicc(mfit.rmb,dat_tps[[i]],3)
			mfit.rbl<-machisplin.aicc(mfit.rbl,dat_tps[[i]],3)
			mfit.mbl<-machisplin.aicc(mfit.mbl,dat_tps[[i]],3)
			mfit.rmbg<-machisplin.aicc(mfit.rmbg,dat_tps[[i]],4)
			mfit.rmbv<-machisplin.aicc(mfit.rmbv,dat_tps[[i]],4)
			mfit.rmbl<-machisplin.aicc(mfit.rmbl,dat_tps[[i]],4)
			mfit.vmbg<-machisplin.aicc(mfit.vmbg,dat_tps[[i]],4)
			mfit.lmbv<-machisplin.aicc(mfit.lmbv,dat_tps[[i]],4)
			mfit.gvlb<-machisplin.aicc(mfit.gvlb,dat_tps[[i]],4)
			mfit.gmlb<-machisplin.aicc(mfit.gmlb,dat_tps[[i]],4)
			mfit.grlb<-machisplin.aicc(mfit.grlb,dat_tps[[i]],4)
			mfit.gvrb<-machisplin.aicc(mfit.gvrb,dat_tps[[i]],4)
			mfit.gvrm<-machisplin.aicc(mfit.gvrm,dat_tps[[i]],4)
			mfit.mvrl<-machisplin.aicc(mfit.mvrl,dat_tps[[i]],4)
			mfit.gmlv<-machisplin.aicc(mfit.gmlv,dat_tps[[i]],4)
			mfit.grlm<-machisplin.aicc(mfit.grlm,dat_tps[[i]],4)
			mfit.grlv<-machisplin.aicc(mfit.grlv,dat_tps[[i]],4)
			mfit.rlbv<-machisplin.aicc(mfit.rlbv,dat_tps[[i]],4)
			mfit.glmrv<-machisplin.aicc(mfit.glmrv,dat_tps[[i]],5)
			mfit.blmrv<-machisplin.aicc(mfit.blmrv,dat_tps[[i]],5)
			mfit.bgmrv<-machisplin.aicc(mfit.bgmrv,dat_tps[[i]],5)
			mfit.bglrv<-machisplin.aicc(mfit.bglrv,dat_tps[[i]],5)
			mfit.bglmv<-machisplin.aicc(mfit.bglmr,dat_tps[[i]],5)
			mfit.bglmr<-machisplin.aicc(mfit.bglmr,dat_tps[[i]],5)
			mfit.bglmrv<-machisplin.aicc(mfit.bglmrv,dat_tps[[i]],6)
			#compare all model fits
			#store all as a vector
			mfit.val<-c(mfit.brt,mfit.rf,mfit.l,mfit.mars,mfit.svm,mfit.gam,mfit.br,mfit.bl,mfit.bm,mfit.lr,mfit.lm,mfit.rm,mfit.bg,mfit.rg,mfit.lg,mfit.mg,mfit.sg,mfit.bv,mfit.rv,mfit.lv,mfit.mv,mfit.rmv,mfit.rbv,mfit.rlv,mfit.rmg,mfit.rbg,mfit.rlg,mfit.rgv,mfit.blv,mfit.bmv,mfit.blg,mfit.bmg,mfit.mlg,mfit.mlv,mfit.mvg,mfit.lvg,mfit.bvg,mfit.rml,mfit.rmb,mfit.rbl,mfit.mbl,mfit.rmbg,mfit.rmbv,mfit.rmbl,mfit.vmbg,mfit.lmbv,mfit.gvlb,mfit.gmlb,mfit.grlb,mfit.gvrb,mfit.gvrm,mfit.mvrl,mfit.gmlv,mfit.grlm,mfit.grlv,mfit.rlbv,mfit.glmrv,mfit.blmrv,mfit.bgmrv,mfit.bglrv,mfit.bglmv,mfit.bglmr,mfit.bglmrv)
			#store names as a vector
			mfit.nam<-c("b","r","l","m","v","g","br","bl","bm","lr","lm","rm","bg","rg","lg","mg","sg","bv","rv","lv","mv","rmv","rbv","rlv","rmg","rbg","rlg","rgv","blv","bmv","blg","bmg","mlg","mlv","mvg","lvg","bvg","rml","rmb","rbl","mbl","rmbg","rmbv","rmbl","vmbg","lmbv","gvlb","gmlb","grlb","gvrb","gvrm","mvrl","gmlv","grlm","grlv","rlbv","glmrv","blmrv","bgmrv","bglrv","bglmv","bglmr","bglmrv")
			mfit<-cbind(as.matrix(mfit.nam),as.matrix(mfit.val))
			colnames(mfit)<-c("id","val")
			sort.mfit<-mfit[order(as.numeric(mfit[,2])),]
			#pick best model
		    f.mod<-sort.mfit[1:1]
			##################################################################################################
			######################################## part 2 run best model(s)  #################################
			##################################################################################################

			#get number of models in best
            ku<-nchar(f.mod)
			#create vector of models
			k=c(1:ku)
		    #split out letters of each model
			mods.run<-unlist(strsplit(f.mod,""))
			pred.elev<- NULL
			res.FINAL<-NULL
			
			for(k in mods.run){	
   			  if (k=="l"){
				mod.run="LM"
				#run model with all points
				mod.lm.tps.FINAL<-stepAIC(lm(mod.form, data = dat_tps[[i]]))
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.lm.tps.FINAL)
				pred.elev<-(pred.elev+pred.elev.2)}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.lm.tps.FINAL)}
				#get residuals
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-dat_tps[[i]]$resp-predict(mod.lm.tps.FINAL)
				res.FINAL<-(res.FINAL+res.FINAL.2)}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-dat_tps[[i]]$resp-predict(mod.lm.tps.FINAL)}
				}
			
  	          ##final models = boosted regression tree	
			  if (k=="b"){
				mod.run="BRT"
				#run model with all points
				mod.brt.tps.FINAL<- gbm.step(data=dat_tps[[i]], gbm.x = 2:(n.covars+1), gbm.y =1, family = "gaussian", tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5, plot.main = FALSE)
				#create raster prediction
				if(is.null(pred.elev)==FALSE){pred.elev.2<- predict(rast_stack, mod.brt.tps.FINAL, n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")
				pred.elev<-(pred.elev+pred.elev.2)}
				if(is.null(pred.elev)==TRUE){pred.elev<- predict(rast_stack, mod.brt.tps.FINAL, n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")}
				#predict at all train sites
				pred.brt.obs <- predict(mod.brt.tps.FINAL, dat_tps[[i]], n.trees=mod.brt.tps.FINAL$gbm.call$best.trees, type="response")
				#get residuals
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-dat_tps[[i]]$resp-pred.brt.obs
				res.FINAL<-(res.FINAL+res.FINAL.2)}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-dat_tps[[i]]$resp-pred.brt.obs}
				}
			
			
  		      ##final models = randomForest	
         	  #if (mfit.rf<mfit.brt & mfit.rf<mfit.lm & mfit.rf=<mfit.mars){
			  if (k=="r"){
				mod.run="RF"
				#run model with all points
				mod.rf.tps.FINAL<- randomForest(mod.form, data = dat_tps[[i]])
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.rf.tps.FINAL, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
				pred.elev<-(pred.elev+pred.elev.2)}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.rf.tps.FINAL, type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)}
				#get residuals
				pred.rf.obs<-predict(mod.rf.tps.FINAL, dat_tps[[i]], type="response", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-dat_tps[[i]]$resp-pred.rf.obs
				res.FINAL<-(res.FINAL+res.FINAL.2)}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-dat_tps[[i]]$resp-pred.rf.obs}
				}
			
   		      ##final models = MARS	
			  if (k=="m"){
				mod.run="MARS"
				#run model with all points
				mod.MARS.tps.FINAL<-earth(mod.form,  data = dat_tps[[i]], nfold=10, ncross=30, varmod.method="lm")
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.MARS.tps.FINAL)
				pred.elev<-(pred.elev+pred.elev.2)}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.MARS.tps.FINAL)}
				#get residuals
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-resid(mod.MARS.tps.FINAL, warn=FALSE)
				res.FINAL<-(res.FINAL+as.vector(res.FINAL.2))}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-as.vector(resid(mod.MARS.tps.FINAL, warn=FALSE))}
				}

   		      ##final models = SVM
			  if (k=="v"){
				mod.run="SVM"
				#run model with all points
				mod.SVM.tps.FINAL<-ksvm(mod.form, data=dat_tps[[i]])
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.SVM.tps.FINAL)
				pred.elev<-(pred.elev+pred.elev.2)}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.SVM.tps.FINAL)}
				#get residuals
				pred.SVM.obs<-predict(mod.SVM.tps.FINAL, dat_tps[[i]])
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-as.vector(dat_tps[[i]]$resp-pred.SVM.obs)
				res.FINAL<-(res.FINAL+res.FINAL.2)}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-as.vector(dat_tps[[i]]$resp-pred.SVM.obs)}
				}

   		      ##final models = GAM
			  if (k=="g"){
				mod.run="GAM"
				#run model with all points
				mod.GAM.tps.FINAL<-gam(mod.form, data=dat_tps[[i]])
				#create raster pred
				if(is.null(pred.elev)==FALSE){pred.elev.2<-predict(rast_stack, mod.GAM.tps.FINAL)
				pred.elev<-(pred.elev+pred.elev.2)}
				if(is.null(pred.elev)==TRUE){pred.elev<-predict(rast_stack, mod.GAM.tps.FINAL)}
				#get residuals
				pred.GAM.obs<-predict(mod.GAM.tps.FINAL, dat_tps[[i]])
				if(is.null(res.FINAL)==FALSE){res.FINAL.2<-as.vector(dat_tps[[i]]$resp-pred.GAM.obs)
				res.FINAL<-(res.FINAL+res.FINAL.2)}
				if(is.null(res.FINAL)==TRUE){res.FINAL<-as.vector(dat_tps[[i]]$resp-pred.GAM.obs)}
				}
			   
			   }
			#wrap up analysis and get final layer
			#divide models by number ensembled
			if(ku>1){
			pred.elev<-(pred.elev/ku)
			res.FINAL<-(res.FINAL/ku)
			}
			
			#calculate Sum of squares and resisdual sum of squares
			rss.m <- sum((res.FINAL) ^ 2)
			tss <- sum((dat_tps[[i]]$resp - mean(dat_tps[[i]]$resp)) ^ 2)
			rsq.model <- 1 - (rss.m/tss) #r squared
			
			##################################################################################################
			###################################### part 3 spline residuals ###################################
			##################################################################################################
			#DETERMINE IF TILING IS NEEDED
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
					}
					#####################################################
					#####################################################
					#feather up
					if (j<nRx){
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
			##################################################################################################
			################################# part 5 finalize results ##################################
			##################################################################################################
			
			#calculate pred at normal scale, suming up kriging pred+res
			pred.elev.i<- brick(pred.elev, final.TPS)
			pred.elev.i.calc<- calc(pred.elev.i, fun=sum)
			
			#extract final values to input points from final raster
            f.actual<-extract(pred.elev.i.calc, dat_tps[[i]][,n.covars:(n.covars+1)])#lat and long input
            
			#calculate Sum of squares and resisdual sum of squares
			rss.final <- sum((dat_tps[[i]]$resp - f.actual) ^ 2)
			rsq.final <- 1 - (rss.final/tss) #r squared
			#out values
			sum.val<-data.frame(out.names[i],f.mod,rsq.model,rsq.final)
			colnames(sum.val)<-c("layer","best model(s):","r2 ensemble:","r2 final:")
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
			return(l)
          }
 
          
#Initiate snowfall
sfInit(parallel=TRUE, cpus=n.cores)
## Export packages
sfLibrary('raster', character.only=TRUE)
sfLibrary('gstat', character.only=TRUE)
sfLibrary('MASS', character.only=TRUE)
sfLibrary('nlme', character.only=TRUE)
sfLibrary('fields', character.only=TRUE)
sfLibrary('dismo', character.only=TRUE)
sfLibrary('randomForest', character.only=TRUE)
sfLibrary('earth', character.only=TRUE)
sfLibrary('kernlab', character.only=TRUE)
sfLibrary('mgcv', character.only=TRUE)

## Export variables
sfExport('dat_tps')
sfExport('rast_stack')
sfExport('i')
sfExport('mod.form')
sfExport('n.covars')
sfExport('out.names')
## Do the run
mySFelevOut <- sfLapply(i, myLapply_elevOut)
## stop snowfall
#sfStop(nostop=FALSE)
sfStop()
f<-mySFelevOut
return(f)
}

##################################################################################################
##################################################################################################
#' Write MACHISPLIN geotiff
#' @param imported table for analysis  
#' @return This tool remove NAs and converts all input data to numeric values for analysis in Humboldt. 
#' @export
#' @examples
#' library(humboldt)
#' 
#' 
#' 
#' # Import a csv as shapefile:
#' Mydata	<-read.delim("sampling.csv", sep=",", h=T)
#' 
#' #load rasters				 
#' ALT = raster("SRTM30m.tif")
#' CURV_V = raster("cs_curv.tif")
#' CURV_H = raster("long_cur.tif")
#' SLOPE = raster("ln_slope.tif")
#' REL_ALT= raster("relative_elevation500m.tif")
#' TWI = raster("TWI.tif")
#' 
#' # function input: raster stack of covarites
#' raster_stack<-stack(ALT,REL_ALT,SLOPE,TWI,CURV_V,CURV_H)
#' 
#' #run an ensemble machine learning thin plate spline 
#' interp.rast<-machisplin.mltps(int.values=Mydata, covar.ras=raster_stack, n.cores=2)
#' 
#' machisplin.write.geotiff(mltps.in=interp.rast)

machisplin.write.geotiff<-function(mltps.in, out.names= NULL){
	n.spln<-mltps.in[[1]]$n.layers
	
	### store rasters
    for (i in 1:n.spln){
		if(i ==1){rast.O<-mltps.in[[i]]$final} else {rast.O<-stack(rast.O, mltps.in[[i]]$final)}
		}
	#write rasters
	if(is.null(out.names)==TRUE){writeRaster(stack(rast.O), names(rast.O), bylayer=TRUE, format='GTiff')}
	if(is.null(out.names)==FALSE){writeRaster(stack(rast.O), names(out.names), bylayer=TRUE, format='GTiff')}

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
	legend3<-"l = linear model; g = generalized additive model (GAM); m = multivariate adaptive regression splines (MARS);"
	legend4<-"v= support vector machines (SVM); r= random forest (RF)"
	legend5<-"NOTE: if R2 Ensemble is greater than R2 Final, then the output model is the ensembled model only (the thin-plate-spline of residuals were not used)"
	write(legend,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend0,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend1,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend2,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend3,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend4,file="MACHISPLIN_results.csv",append=TRUE)
	write(legend5,file="MACHISPLIN_results.csv",append=TRUE)
}

