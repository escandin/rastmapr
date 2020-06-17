namaskplus <- function(x, toval, otherval) {
  
  # convertes NAs to to a value indicated by toval and any other values
  # to a new value indicated by otherval
  require(raster)
  v <- raster::getValues(x)
  isna=is.na(v)
  isnotna=!is.na(v)
  v[isna]=toval
  v[isnotna]=otherval
  x <- raster::setValues(x, v)
  return(x)
}
namask <- function(x) {
  ###### DEPRICATE!
  # convertes NAs to 1s and no NAs to NAs
  require(raster)
  v <- raster::getValues(x)
  v[v>0]=0
  v[is.na(v)]=1
  v[v==0]=NA
  x <- raster::setValues(x, v)
  return(x)
}
maskfun <- function(x, thresh, aboveval, belowval) {
  # sets a threshold value (thresh). Any pixels above 
  # thresh are converted to aboveval
  # Any pixels bellow or equal to thresh are converted to belowval
  require(raster)
  v <- raster::getValues(x)
  v[v>thresh]=aboveval
  v[v<=thresh]=belowval
  x <- raster::setValues(x, v)
  return(x)
}
stackmask= function(inrast=inrast, maskrast=TRUE){
  # THIS FUNCTION HAS TO BE OPTIMIZED. IT TAKES TOO LONG
  # masks out any pixels that have NAs in at least 1 band 
  # in a raster stack
  # maskrast: TRUE applies the mask to the raster stack
  # output: If maskrast=true a list with the stack masked and
  # the mask produced. If maskrat=FALSE, the raster masck)
  msk<-max(inrast)
  msk<-maskfun(msk, 0,1,NA)
  #inrast<-raster::mask(inrast, msk)
  if (maskrast==TRUE){
    inrast<-raster::mask(inrast, msk) # this is more than a minute faster than using the raster::mask function
    out=list(inrast,msk)} else {
      out=msk}
  return(out)
}
commonExtent=function(inlist=rasterlist, method="intersection"){
  # calculates the common extent between all the raster objects within a list.
  # The objects should all overlap at least partially and have the same projection, spatial resolution and be aligned
  # This can be written in a function. It can calculate either the largest extent that includes
  # all the listed elements (union) orthe minimum common extent (intersection).
  for (i in 1: length(inlist)){
    ext=as.matrix(extent(inlist[[i]]))
    if(i==1){
      extall=ext} else {extall=rbind(extall, ext)}
  }
  
  extall=data.frame(extall)
  if (method=="intersection"){
  minx=max(extall$min[c(seq(1,nrow(extall), by=2))])
  maxx=min(extall$max[c(seq(1,nrow(extall), by=2))])
  miny=max(extall$min[c(seq(2,nrow(extall), by=2))])
  maxy=min(extall$max[c(seq(2,nrow(extall), by=2))]) 
  } else if (method=="union"){
    minx=min(extall$min[c(seq(1,nrow(extall), by=2))])
    maxx=max(extall$max[c(seq(1,nrow(extall), by=2))])
    miny=min(extall$min[c(seq(2,nrow(extall), by=2))])
    maxy=max(extall$max[c(seq(2,nrow(extall), by=2))])
  } else {stop("unvalid method name")
    }
  e=extent(minx, maxx, miny, maxy)
  return(e)
}
crossvalidsp=function(dataset=indata, classcolname="classes", validprop=0.2){
  require(maptools, raster)
  uniqueclass=unique(dataset[[classcolname]])
  for (i in 1:length(uniqueclass)){
    classrows=which(dataset[[classcolname]]==uniqueclass[i],arr.ind=TRUE)
    classpoly=dataset[classrows,]
    validrows=sample(c(1:nrow(classpoly)), round(nrow(classpoly)*validprop))
    validpoly=classpoly[validrows,]
    calibpoly=classpoly[-validrows,]
    if (i==1){
      validtot=validpoly
      calibtot=calibpoly} 
    else {
      validtot=spRbind(validtot, validpoly)
      calibtot=spRbind(calibtot, calibpoly)
    }
  }
  return(list(calibtot, validtot))
}
conmatrixmerge=function(reclassx=conmatrix, rcl=fromtomatrix){
  # This function merges pairs of classes in an accuracy matrix
  # and recalculates accuracy metrics.
  # reclassx= the input accuracy matrix as a table
  # rcl= a nx2 reclassification matrix with a number representing 
  # the target class object in the first column and the class to be merged to
  # in the second column
  for (i in 1:nrow(rcl)){
    reclassx[rcl[i,2],]=reclassx[rcl[i,2],]+reclassx[rcl[i,1],]
    reclassx[,rcl[i,2]]=reclassx[,rcl[i,2]]+reclassx[,rcl[i,1]]
  }
  reclassx=reclassx[,-rcl[,1]]
  reclassx=reclassx[-rcl[,1],]
  overall=sum(diag(reclassx))/sum(reclassx)
  producers=diag(reclassx)/apply(reclassx, 2,sum)
  users=diag(reclassx)/apply(reclassx, 1,sum)
  accustats=cbind(users,producers)
  out=list(overall, accustats, reclassx)
  names(out)=c("overall","users/producers", "matrix")
  return(out)
}
AccuStats=function(prefix="class2013_AP_raw_iter", suffix="accumatrix.RData", iter=30){
  # extracting accuracy results from confusion matrices. Add this to accuplot!
  for (i in 1:iter){
    outname=paste(paste(prefix, i, sep="_"), suffix, sep="_") 
    load(outname)
    accuresults=data.frame(conmatrix$byClass, stringsAsFactors=FALSE)
    users=accuresults$Precision
    producers=accuresults$Recall
    overall=conmatrix$overall[[1]]
    classnames=rownames(accuresults)
    classnames=str_remove(classnames, "Class: ")
    if (i==1){
      allusers=users
      allproducers=producers
      alloverall=overall} else {
        allusers=data.frame(cbind(allusers, users), stringsAsFactors=FALSE)
        allproducers=data.frame(cbind(allproducers, producers), stringsAsFactors=FALSE)
        alloverall=rbind(alloverall, overall)
      }
  }
  allusers=data.frame(cbind(classnames, allusers), stringsAsFactors=FALSE)
  allproducers=data.frame(cbind(classnames, allproducers), stringsAsFactors=FALSE)
  alloverall=alloverall
  mean(alloverall)
  sd(alloverall, na.rm=TRUE)
  allusers=data.frame(allusers, stringsAsFactors=FALSE)
  allproducers=data.frame(allproducers,stringsAsFactors=FALSE)
  
  # summarize user's and producer's mean and sd accuracy in a table
  finalstats=data.frame(cbind(round(apply(allusers[,2:ncol(allusers)], 1,  mean), digits=2), 
                              round(apply(allusers[,2:ncol(allusers)], 1,  sd), digits=2), 
                              round(apply(allproducers[,2:ncol(allproducers)], 1,  mean), digits=2),
                              round(apply(allproducers[,2:ncol(allproducers)], 1,  sd), digits=2)), stringsAsFactors=FALSE)
  # I have to do this twice otherwise data.frame converts numbers into factors.
  finalstats=cbind(classnames,finalstats)
  names(finalstats)=c("classnames", "users mean", "users sd", "producers mean", "producers sd")
  
  return(list(allusers, allproducers, alloverall, finalstats))
}
stack2df=function(inrast=rastack,invec=polygons, classcolname="class"){
  # extracts into a data frame the pixel values for all bands from different classes
  # defined in a spatial dataframe
  # inrast: the raster dataset containing pixel values to extract in [[1]]
  # invec: spatial dataframe object defining the locations where the data
  # should be extracted from 
  # classcolname: the column in the spatial dataframe containing the names 
  # of the attributes associated to those areas
  # value: a data frame with columns representing the pixel values in each band for
  # the areas labeled as defined by classcolname
  if (is.null(intersect(extent(invec), extent(inrast)))){
    stop("the extents of inrast and invec do not overlap")
  }
  if(as.character(crs(inrast))!= as.character(crs(invec))){
    stop("inrast and invec should have the same projection")
  }
  # required function
  extractval=function(inraster=inrast, msk=msk){
    outvector=raster::mask(inraster, msk) 
    outvector=na.omit(raster::getValues(outvector))
    return(outvector)
  }
  
  # assign class ID to each class
  invec$class_ID=rep(NA, nrow(invec@data))
  for (i in 1:length(invec[[classcolname]])){
    invec$class_ID[which(invec[[classcolname]]==levels(invec[[classcolname]])[i])]=i
  }
  
  # mask the input raster including  pixels with valid values in all bands only
  inrast=stackmask(inrast)
  
  # create a raster of class_ids. TRY gdalUtils::gdal_rasterize. It might be faster!!!
  calibrast=raster::rasterize(invec, inrast[[2]], field=invec$class_ID)
  calibmsk<-maskfun(calibrast, 0, 1, NA)
  calibmsk=raster::mask(calibmsk, inrast[[2]])
  
  # Extract pixel values into a dataframe
  class_ID=(extractval(calibrast, calibmsk))
  dataset=data.frame(matrix(data=NA, nrow=length(class_ID), ncol=nlayers(inrast[[1]])))
  
  # add a column with a class name
  dbclassname=rep(NA, length(class_ID))
  for (i in 1:length(levels(invec[[classcolname]]))){
    dbclassname[which(class_ID==i, arr.ind=TRUE)] = levels(invec[[classcolname]])[i]
  }
  commonclasses= match(sort(unique(dbclassname)), sort(levels(invec[[classcolname]])))
  if(length(commonclasses)< length(levels(invec[[classcolname]]))){
    missing=sort(levels(invec[[classcolname]]))[-commonclasses]
    warning(paste(paste("the class", missing, sep= " "), 
                  "has no valid pixels in input raster", sep=" "))
    print(paste(paste("Warning: the class", missing, sep= " "), 
                "has no valid pixel values in the input raster", sep=" "))
  }
  
  dataset=cbind(class_ID, dbclassname, dataset)
  rm(class_ID, dbclassname)
  
  #I HAVE TO REMOVE DATA THAT HAS BEEN MASKED FROM THE RASTER
  # I MIGHT HAVE TO BRING BACK THE STACKMASK FUNCTION
  # ranges=c(-Inf,1,NA, 1,99999,1)
  # calibmsk=reclassify(calibrast,ranges)
  # calibmsk=calibmsk*inrastmskd[[2]]
  
  for (i in 1:nlayers(inrast[[1]])){
    dataset[,i+2]=extractval(inrast[[1]][[i]], calibmsk)
    print(paste(i, "layers extracted", sep=" "))
  }
  names(dataset)=c("class_ID", "class_name", names(inrast[[1]]))
  return(dataset)
}
AccuPlot=function(conmatrix=confusion, location="top"){
  #Plots users and producers accuracy based on the results of a confusion matrix 
  # obtained through the application of the function caret::confusionMatrix
  # location: location of the legend. Options:"bottomright", "bottom", "bottomleft",
  # "left", "topleft", "top", "topright", "right"  "center"
  require(stringr)
  accuresults=data.frame(conmatrix$byClass, stringsAsFactors=FALSE)
  users=round(accuresults$Precision, digits=3)
  producers=round(accuresults$Recall, digits=3)
  overall=conmatrix$overall[[1]]
  classnames=rownames(accuresults)
  classnames=str_remove(classnames, "Class: ")
  #names(accuresults)=classnames
  par(mar = c(7, 4, 2, 2) + 0.2)
  barplot(rbind(users,producers),col=c("aquamarine3","coral"), 
          names.arg=classnames,  beside = TRUE, ylab= "accuracy (%)", las=2, cex.names=0.7)
  legend(location, legend=c("Users", "Producers"), box.lty=0, bg= "transparent",
         col=c("aquamarine3","coral"), lty=1:2, cex=0.8)
  # as.numeric(levels(aRnStats[[4]][,2]))
  #
  #plot(accuplot)
  classaccuracy=data.frame(cbind(as.numeric(users), as.numeric(producers)), 
                           stringsAsFactors=FALSE)
  # I have to do this twice otherwise data.frame converts numbers into factors.
  classaccuracy=cbind(classnames, classaccuracy, stringsAsFactors=FALSE)
  names(classaccuracy)=c("class name", "users", "producers")
  return(list(overall, classaccuracy))
}
rasterCCA=function(xstack, ystack, propsamp=1){
  # applies a canonical correlation transformation to two raster stacks
  # THIS FUNCTION NEEDS TO BE MORE EFFICIENT. I have to figure out how
  # to do CCA directly in the raster and ignore NAs.
  # SELECTING A PROPSAMP DIFFERENT THAN 1 PRODUCES SUBSTANTIALLY DIFFERNT RESULTS
  if(nlayers(xstack)!=nlayers(ystack)){
    stop("the xstack and ystack stacks should have the same number of bands")
  }
  msk=mask(stackmask(xstack, maskrast=FALSE),
           stackmask(ystack, maskrast=FALSE))
  xstack=mask(xstack, msk)
  ystack=mask(ystack, msk)
  
  xm=getValues(xstack)
  ym=getValues(ystack)
  if(is.numeric(propsamp) & propsamp>=0 & propsamp<1){
    xtrain=raster::sampleRandom(xstack, size = nrow(xstack)*ncol(xstack)*propsamp, na.rm = TRUE)
    ytrain=raster::sampleRandom(ystack, size = nrow(ystack)*ncol(ystack)*propsamp, na.rm = TRUE)
  } else if (propsamp==1){
    xtrain=na.omit(xm)
    ytrain=na.omit(ym)
  } else {
    stop("propsamp should be a number between 0 and 1")
  }
  CCA=cancor(xtrain,ytrain)
  xCCA=xm %*% CCA$xcoef
  yCCA=ym %*% CCA$xcoef
  xstack=setValues(xstack, xCCA)
  ystack=setValues(ystack, yCCA)
  return(list(CCA, xstack, ystack))
}
smg=function(inlist=listr, cca="FALSE", method="none", refitem=NA, mosaicitems=c(1,3), 
             QAbandname="pixel_qa", cloudbuff=NA, savelcloudmsk=TRUE,
             indem=NA, normbands=seq(1, 7), pval.pif=1e-02, pval.chg=0.99,
             sensor="OLI", verbose=TRUE, norm.ext=NULL){
  # smg (seamless mosaic generator) produces a cloud masked and radiometrically normalized 
  # raster mosaic from landsat images.
  # inlist: list of input rasters
  # method: "none" for no normalization, "aRn" for std diff normalization or "cor" for pearson correlation (pending)
  # refitem= item in the list corresponding to the raster to use as reference 
  # for the radiometric normalization. This image should overlap at least partially
  # with all mosaicitems in order to work.
  # mosaicitems: items to be normalized and mosaiced. If mosaicitems=1, no mosaics are built
  # QAband (inherited): position of the landsat QA band in each rasterstack contained in 
  #         inlist. If QA=NA, no cloud masking is performed
  # sensor (inherited): Sensor to encode. Options: c("OLI", "TIRS", "ETM+", "TM", "MSS")
  # transform: whether to apply a canonical transformation to the data before normalizing
  # required functions and packages
  require(rgdal, raster, RStoobox, aRn)
  is.integer0 <- function(x){
    is.integer(x) && length(x) == 0L
  }
  validmethods=c("s3d", "ed", "cor", "sam", "none")
  if (is.integer0(which(validmethods==method))){
    stop("please provide a valid method: s3d, ed, cor, sam, none")
  }
  
  stacknames=names(inlist)
  maskitems<-unique(c(refitem, mosaicitems))
  maskitems<-maskitems[!is.na(maskitems)]
 
  # perform cloud masking
  if(!is.na(QAbandname)){
    # select items to mask from the list
    #bnames=names(inlist[[refitem]])
    #cloudmasking:
    print('masking clouds')
    for(i in maskitems){
      layernames=names(inlist[[i]])
      msk=RStoolbox::classifyQA(inlist[[i]][[QAbandname]], sensor=sensor) #do Sys.time to know whether it's faster with get/setValues
      if(is.numeric(cloudbuff)){
         msk=buffer(msk, width=cloudbuff)}
      msk=namaskplus(msk, 1, NA)
      #msk=namask(msk) # see function above
      if(savelcloudmsk==TRUE){
      writeRaster(msk, filename=paste(paste("CloudMsk", stacknames[i], sep="_"), "tif", sep="."))
      }
      inlist[[i]]=inlist[[i]][[normbands]] # only select bands to normalize
      inlist[[i]]=raster::mask(inlist[[i]], msk) ## this is more than a minute faster than using the raster::mask function
      names(inlist[[i]])=layernames[normbands]
       #names(inlist[[i]])=bnames[normbands]
      print(paste(paste("item", i, sep=" "), "masked", sep=" "))
    }
    rm(msk)
  }
  
  # perform topographic correction
  if(class(indem)[1]=="RasterLayer"){
    print ("topographic correction started")
    for (j in maskitems){ 
      # ensure that the dem is aligned with the raster stack
      if(as.character(crs(indem))!=as.character(crs(inlist[[j]]))){
        print("reprojecting raster dem")
        indem<-raster::projectRaster(indem, inlist[[j]], 
                                     datatype=dataType(inlist[[j]][[normbands[1]]]))
      } else if(as.character(extent(indem))!=as.character(extent(inlist[[j]]))){
        print("resampling raster dem") 
        # OJO: dataType is not working. It only allows 
        subdem<-raster::resample(indem, inlist[[j]], 
                                 datatype=dataType(inlist[[j]][[normbands[1]]]))
      } else {
        subdem<-indem
      }
      # This is a workaround to change the data type. 
      # The format of the file changes to FLT4S after resampling regardless of
      # whether I add the datatype argument or change the global rasterOptions 
      # writeRaster(subdem, paste(tmpDir(), "subdem.tif", sep=""), 
      #            datatype=dataType(inlist[[j]][[normbands[1]]]))
      #subdem=raster(paste(tempdir, "subdem.tif", sep="/"))
      print("producing common mask between dem and raster stack")
      msk<-stackmask(stack(inlist[[j]], subdem), maskrast=TRUE)
      inlist[[j]]<-msk[[1]][[1:nlayers(stack(inlist[[j]]))]]
      subdem<-msk[[1]][[nlayers(stack(inlist[[j]]))+1]]
      rm(msk)
      #print("extracting topographic variables")
      #slope<-raster::terrain(dem, opt='slope')
      #aspect<-raster::terrain(dem, opt='aspect')
      #terrainvar<-raster::stack(slope,aspect)
      #names(terrainvar)<-c('slope', 'aspect')
      #rm(slope,aspect,dem)
      meta<-RStoolbox::readMeta(paste(names(inlist)[[j]], "xml", sep="."))
      print ("performing topographic correction")
      inlist[[j]]<-RStoolbox::topCor(inlist[[j]], dem=subdem, metaData=meta, method="C")
      rm(subdem)
      # THIS PART OF THE WORKAROUND: remove temporary dem file so that a new one can be saved
      #file.remove(paste(tempdir, "subdem.tif", sep="/"))
      print(paste(paste("item", j, sep=" "), "topo-corrected", sep=" "))
    }
  }

  # normalize
  normitems=mosaicitems[which(mosaicitems!=refitem)]
  if(method=="s3d"){
    print("normalizing")

    for (i in normitems){
      PIFlist=c(inlist[[i]],inlist[[refitem]])
      names(PIFlist)=names(inlist)[c(i, refitem)]
      models <- s3d(strips=PIFlist, cca=cca)
      # make sure that the items in the list have a name
      save(models, file=paste(paste("s3dModels", stacknames[i], sep="_"), "RData", sep=".")) 
      #save(models, file=paste("modelsItem", i, sep="_"))
      if(verbose == TRUE) {print(summary(models))}
      names(models$parameters)=names(PIFlist[[1]])
      inlist[[i]]<-aRn::CalibrateRaster(models$parameters, PIFlist)
      # if(verbose == TRUE) {plot(models)} # I stoped this to reduce processing time
      print(paste(paste("item", i, sep=" "), "normalized", sep=" "))
    }
  } else if(is.integer0(which(validmethods[2:4]==method))==FALSE){
    print("normalizing")
    for (i in normitems){
      normobj<-RStoolbox::pifMatch(inlist[[i]], inlist[[refitem]], method=method, 
                                   quantile=1-pval.pif, returnModels=TRUE)
      save(normobj, file=paste(paste(paste(method, "Models", sep=""), stacknames[i], sep="_"), "RData", sep=".")) 
      inlist[[i]]<-normobj[[1]]
      #if(verbose == TRUE) {plot(models)}
      print(paste(paste("item", i, sep=" "), "normalized", sep=" "))
    }
  }
  
  # mosaic TRY THE FUNCTION gdalUtils::mosaic_rasters IT MIGHT BE FASTER!!
  mosaicYr=inlist[[mosaicitems[1]]]
  if (length(mosaicitems)>1){
    print("mosaicking layers")
    mosaicYr=raster::merge(inlist[[mosaicitems[[1]]]], inlist[[mosaicitems[[2]]]])
    print("2 layers mosaicked")
  }
  
  if(length(mosaicitems)>2){
    for (i in 3:length(mosaicitems)){
      mosciacYr=raster::merge(mosaicYr, inlist[[mosaicitems[[i]]]])
      print(paste(i, "layers mosaicked", sep=" "))
    }
  }
  names(mosaicYr)=names(inlist[[refitem]])[normbands]
  return(mosaicYr)
}
tmg=function(inrast=x, calibObject=Calib, ntrees=1000, classcolname="class", 
             validObject=Valid, maptype="class", savefiles=TRUE, 
             outfile="outfilename", plot=TRUE, tmproduce=TRUE){
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
  require(raster, randomForest, rgdal, caret)
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
    if(length(levels(validObject[[classcolname]]))!=length(rownames(class.RF$confusion))){
      stop('The number of classes in validObject should be the same as the number of classes in calibObject')
    }
    print("extracting pixel values for validation areas")
    if(is.na((mean(match(rownames(class.RF$confusion),levels(validObject[[classcolname]])))))==TRUE){
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
      writeRaster(tm, filename=paste(outfile, "tif", sep="."))
      }
  }
  # return(tm) # not sure what object to return since it depends on the inputs
}
lcupdate=function(refclass=NA, tarclass=NA, nochmsk=NA){
  # lcupdate produces an update of a target land cover map (tarclass) by replacing
  # the classification values of pixels located in areas considered as no change as
  # denoted in a no change mask (nochmsk). See s3d for methods to produce no change masks).
  # 
  # refclass: reference classification
  # tarclass: target classification
  # nochmsk: #mask representing areas labeled as no change between refclass and tarclass
  chmsk=namaskplus(msk, 1, NA)
  #chmsk=namask(nochmsk) 
  refclass=raster::mask(refclass, nochmsk)
  tarclass=raster::mask(tarclass, chmsk)
  tarclass=raster::merge(tarclass,refclass)
  return(tarclass)
}


