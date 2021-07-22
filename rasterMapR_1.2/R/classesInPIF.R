classesInPIF <- structure(function #Classes in PIF
### This function is under testing and has not been implemented
###yet. For computing \code{PIF} see \code{\link{PIFmodel}},
###\code{\link{thresraster}}, \code{\link{nochg}}, and
###\code{\link{calibrationParameters}}. The function calculates
###classes into brackets of Pseudo-Invariant Features.
                          ##details<< Selected \code{PIF} are
                          ##scaled by computing the square of the
                          ##scaled differences between two rasters,
                          ##see \code{\link{scale}}. The scaled
                          ##\code{PIF} are then averaged using classes
                          ##of reflectance in the compared images
                          ##(\code{PIF} classes). The classes are
                          ##computed using \code{\link{quantile}}.
(
    pifdt, ##<<\code{data.frame}. Data set of PIF
    brackets = 15 ##<<\code{numeric}. Number of brackets.
) {
    if(is.null(pifdt))
        return(NULL)
    pifdtna <- na.omit(pifdt)
    pifdtna <- pifdtna[order(pifdtna$'x'),]
    x2ranges <- function(x, brackets){
        rang <- range(x, na.rm = TRUE)
        rn. <- do.call('seq', as.list(rang))
        pr. = seq(0, 1, by = 1/brackets)
        brk <- stats::quantile(rn., probs = pr.)
        brk <- as.numeric(brk)
        
        ## brlngt <- diff(rang)/n
        ## rn <- 1:(n - 1)
        ## sqrn <- rn * brlngt + rang[1L] 
        ## brkrng <- c(rang, sqrn)
        ## brk <- brkrng[order(brkrng)]
        
        ## Note: last chunck of code is based on
        ## nbrackets <- n
        ## indata <- pifdtna
        ## xval <- pifdtna$'x'
        ## yval <- pifdtna$'y'
        ## bracketlength=(max(xval)-min(xval))/nbrackets
        ## bracketbins=seq(1,nbrackets-1)
        ## bracketranges=c(min(xval),
        ## bracketbins*bracketlength+min(xval), max(xval))
        ## brk <- bracketranges
        
        return(brk)}
    range. <- x2ranges(pifdtna$'x', brackets)
    cuts <- cut(pifdtna$'x', breaks = range.)
    dfchunk <- split(pifdtna, cuts)
    ## row.names(pifdtna) <- 1:nrow(pifdtna) 
    ## dfchunk <- split(pifdtna, factor(sort(rank(row.names(pifdtna))%%n)))
    apd <- lapply(dfchunk, function(x)
        apply(x[,c('y','x')], 2, mean))
    asd <- lapply(dfchunk, function(x)
        sd(x$'y'))
    apsd <- Map(function(x,z)c(x, sd = z),
                apd, asd)
    apd <- data.frame(na.omit(do.call('rbind', apsd)))
    return(apd)
### \code{data.frame}.
} , ex=function(){
    classesInPIF(NULL)
})
