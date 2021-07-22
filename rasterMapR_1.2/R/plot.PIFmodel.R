plot.PIFmodel <- structure(function #PIF-model plots
### This function plots linear models of Pseudo-Invariant Features
### (\code{PIF} models).
(
    x, ##<<\code{list}. Set of \code{PIF} models such as that produced
       ##by \code{\link{PIFmodel}}.  
    ... ##<< Additional arguments in \code{\link{plot}} other than
        ##\code{axes}, \code{col}, \code{ylab}, \code{xlab}, or
        ##\code{main}
){
    rmse <- function(m, mm){
        preds <- predict.lm(m, mm)
        sqd <- (mm$'y' - preds)^2
        smres <- sqrt(mean(sqd, na.rm = TRUE))
        smres <- round(smres, 3)
        sqtxt <- paste('RMSE =', smres)
        return(sqtxt)}

    raw. <- as.list(attr(x, 'env'))
    rawLs <- raw.[['raw']]
    labrm <- Map(function(l, u)
        rmse(l, u), x, rawLs)

    par(mfrow = n2mfrow(length(rawLs)))
    par(oma = c(4, 4, 1, 1))
    par(mar = c(1.5, 1.5, 3, 1.5))
    for(i in 1:length(x)){
        plot(y ~ x, rawLs[[i]], axes = FALSE,
             col = '#00000033',
             ylab = '', xlab = '',
             main = names(rawLs[i]), ...)
        box('plot', col = 'grey50', lwd = 2)
        axis(1, col = 'grey50',lwd = 2)
        axis(2, col = 'grey50',lwd = 2)
        abline(x[[i]], col = 'red')
        points(y ~ x, x[[i]]$'model',
               col = 'red', cex = 0.6, pch = 19)
        usr <- par('usr')
        text(usr[1], usr[4], labrm[[i]],
             adj = c(-0.25,3), cex = 0.8)
    }
    
    nms. <- raw.[['nm']]
    fch <- function(nm.., fr = 0.4){
        nch <- round(fr * nchar(nm..), 0)
        pch <- substr(nm.., 1, nch)
        return(paste(pch, '...'))}
    mtext(fch(nms.[2]), side = 1,
          outer = TRUE, font = 2,
          line = 1, col = 'grey50')
    mtext(fch(nms.[1]), side = 2,
          outer = TRUE, font = 2,
          line = 1, col = 'grey50')
    par(oma=c(1.5,0,0,0))
    shortname <- "aRn package" # or maybe a filename
    mtext(paste(shortname, " ",
                format(Sys.time(), "%Y-%m-%d %H:%M")),
          cex=0.7, line=0, side=SOUTH<-1,
          adj=0.05, outer=TRUE, col = '#00000033')
    par(mfrow=c(1,1))
### \code{list}. Set of \code{\link{lm}} models
} , ex=function(){
        tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
    strips <- RasterIntersection(stack)
    ## thrs <- thresraster(strips[[2L]], strips[[1L]])
    ## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
    ## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
    model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
    ## plot(model)
})
