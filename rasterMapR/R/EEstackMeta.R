EEstackMeta <- structure(function #Stack EE data
### This function is under development. The function implements
### \code{\link{stackMeta}} to stack Landsat hyper-spectal layers
### using embedded metadata. Use \code{\link{EEstackWithoutMeta}}
### instead.
(
    tarPath, ##<<\code{character} or \code{List}. Name(s) of the
             ##commpresed data set(s). These can include the file
             ##paths.
    bands = NULL, ##<<\code{numeric}. Bands to be used.
    ... ##<< Other arguments in \code{\link{stackMeta}}
){
    fn2 <- function(tarPath, ...){ 
        if(is.null(tarPath))
            return(NULL)
        tarNm <- basename(tarPath)
        Nm <-  sub('.tar|.tar.gz','',tarNm)
        temp. <- tempdir()
        untar(tarPath, exdir = temp.)
        dbase <- list.files(temp.)
        pathtifs <- file.path(temp., dbase)
        tif <- pathtifs[grepl('MTL.txt', pathtifs)]
        tifs <- tif[grepl(Nm, tif)]
        meta <- readMeta(tifs)
        stm <- tryCatch(stackMeta(meta, ...),
        error = function(e) NULL)}
    ## if(length(tarPath) == 1){
    ##     mta <- fn2(tarPath, ...)
    ## }else{
        mta <- Map(function(x, ...)
            fn2(x, ...), tarPath, ...)
## }
    ## if(!is.null(bands))
    ##     mta <- lapply(mta, function(x)
    ##         x[[bands]])
    return(mta)
### RasterBrick.
} , ex=function(){
    EEstackMeta(NULL)
})
