tmg <- structure(function #Thematic map generator
### This function produces a thematic map based on either training
### data or previously calibrated randomForest object.
                       ##details<< 
(
    inrast, ##<<\code{list}. List of input rasters
    calibObject, ##<<\code{logical}. 
    ntrees=1000, ##<<\code{numeric}. "none" for no normalization,
                   ##"aRn" for std diff normalization or "cor" for
                   ##pearson correlation (pending)
    classcolname="class", ##<<\code{}. Item in the list corresponding
                          ##to the raster to use as reference
    validObject, ##<<\code{numeric}. Items to be normalized and
                        ##mosaiced. If mosaicitems=1, no mosaics are
                        ##built.
    maptype="class", ##<<\code{character}
    savefiles=TRUE, ##<<\code{}.
    outfile="outfilename", ##<<\code{character}.
    plot=TRUE, ##<<\code{logical}.
    tmproduce=TRUE, ##<<\code{numeric}.
    dt="INT1U" ##<<\code{character}.
){
  # tmg (thematic map generator) produces a thematic map based on either training data or
  # previously calibrated randomForest object. 
  # inrast: input raster with each band corresponding to a predictor
  # calibObj: either a spatial dataframe containing the training data or a previously
  # calibrated randomForest object
  # classcolname: the column in the attributes of the spatial dataframe containing the name
  # of the classes to predict (OJO: generalize for random forest regression)
  # validdata: a spatial dataframe containing the validation data. If validdata is no spatial dataframe, no validation is performed
  # maptype: whether a randomForest regression or classification will be performed
  # (OJO look for the name to be used when it is a random forest regression. I believe it is "response") )
  ## require(raster, randomForest, rgdal, caret)
  names(inrast)=paste("b", c(1:length(names(inrast))), sep="")
  
  # if calibobject is a spatial object (point, polygon, line), then use it to calibrate
  # a random forest classification
  if((substr(class(calibObject)[1],1,7)=="Spatial")==TRUE){
    #if((class(calibObject)=="data.frame")==TRUE){
    
    # Make the band names and the names of the variables in the RF object compatible
    # OJO: These names are reverted to the original band names when I use stack2df()
    
    # extract pixel values in a dataframe
    print("extracting pixel values for training areas")
    calibdf=stack2df(inrast=inrast, invec=calibObject, classcolname=classcolname)
    #calibdf=calibObject
    # for some reason satack2df is not assigning the names of inrast to the respective columns in calidbdf
    names(calibdf)[3:ncol(calibdf)]<- names(inrast) 
    calibdf[[2]]=as.factor(calibdf[[2]])
    
    # produce a formula with the band names to be entered in the RF classification equal
    # to the input band names
    formulaStr=as.formula(paste(paste(names(calibdf)[2],"~"), paste(names(inrast), 
                                                                    sep=" ", collapse="+"), sep=" "))
    # produce classification object
    print("calibrating classification algorithm")
    class.RF<-randomForest::randomForest(formula = formulaStr, data=calibdf, importance = TRUE, 
                                         ntree =ntrees, proximity=FALSE)
    if(savefiles==TRUE){
      save(class.RF, file=paste(outfile, "_RF.RData", sep=""))
    }
  } else if(substr(class(calibObject)[1],1,12)=="randomForest"){
    # if the object is randomForest type, then us it to classify inraster
    class.RF<-calibObject
    
    #test that the number of predictors in the RF object is 
    # equal to the number of bands in inrast
    # OJO: We have to make sure that the predictors in the RF object 
    # should have the same names as names(inrast)
    if(length(names(inrast))!=length(attributes(class.RF$terms)$term.labels)){
      stop("the number of bands in the input raster should be the same as the number of predictors to calibrate the random forest object")
    }
  } else {
    stop("calibObject is the wrong class. Provide either a spatial dataframe or randomForest object")        
  }
  
  if((substr(class(validObject)[1],1,7)=="Spatial")==TRUE){
    if(length(unique(validObject[[classcolname]]))!=length(rownames(class.RF$confusion))){
      stop('The number of classes in validObject should be the same as the number of classes in calibObject')
    }
    print("extracting pixel values for validation areas")
    if(is.na((mean(match(rownames(class.RF$confusion),unique(validObject[[classcolname]])))))==TRUE){
      stop("The names of classes in the calibration object do not match the ones in the validation object")
    }
    
    validdf=stack2df(inrast=inrast, invec=validObject, classcolname=classcolname)
    #validdf<-validObject
    # for some reason satack2df is not assigning the names of inrast to the respective columns in validdf
    names(validdf)[3:ncol(validdf)]<- names(inrast)
    
    print("generating predictions for validation areas")
    results<-raster::predict(class.RF, validdf, type=maptype)
    print("calculating accuracy metrics")
    conmatrix<- caret::confusionMatrix(as.factor(results),as.factor(validdf$class_name))
    if(savefiles==TRUE){save(conmatrix, file=paste(outfile, "accumatrix.RData", sep="_"))}
  }
  
  if (tmproduce==TRUE){
    print("classifying map")
    tm<-raster::predict(object=inrast, type=maptype, model=class.RF, 
                        fun=predict, ext=NULL, const=NULL, na.rm=TRUE)
    if(savefiles==TRUE){
      print("saving classified map to disc")
      writeRaster(tm, filename=paste(outfile, "tif", sep="."), datatype=dt)
      }
  }
  ### return(tm) # not sure what object to return since it depends on the inputs
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1:4))
    strips <- RasterIntersection(stack)
    ## thrs <- thresraster(strips[[2L]], strips[[1L]])
    ## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
    ## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
    ## model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    ## calib <- CalibrateRaster(model, stack)
    ## merged <- merge(calib, stack[[2L]][[names(calib)]])
    ## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
})
