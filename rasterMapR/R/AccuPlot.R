AccuPlot <- structure(function #Plot users and producers accuracy 
                      ### This function plots users and producers accuracy based on the
                      ### results of a confusion matrix obtained through the application of
                      ### the function caret::confusionMatrix.
                      ##details<< Location: location of the
                      ##legend. Options:"bottomright", "bottom",
                      ##"bottomleft","left", "topleft", "top",
                      ##"topright", "right" "center"
                      (
                        conmatrix, ##<<\code{}. 
                        location="top" ##<<\code{}. 
                      ){
                        ## require(stringr)
                        colsum=colSums(supervised$validation$performance$table)
                        rowsum=rowSums(supervised$validation$performance$table)
                        diagonal=diag(supervised$validation$performance$table)
                        producers=round(diagonal/colsum, digits=3)
                        users=round(diagonal/rowsum, digits=3)
                        overall=sum(diagonal)/sum(colsum)
                        classnames=names(users)
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
                        ### \code{}... 
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
