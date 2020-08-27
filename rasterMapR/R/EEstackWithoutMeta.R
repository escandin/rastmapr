EEstackWithoutMeta <- structure(function #Stack EE avoiding metadata
### This function stacks Landsat data shared by Earth Explorer while
### avoiding embedded metadata.
                                ##details<< Compressed files of Multi-
                                ##or hyper-spectral layers
                                ##(\code{.tar}) are recursively
                                ##decompressed, stacked, and
                                ##renamed. The data is processed using
                                ##the names of the layers while
                                ##ignoring the metadata. However, a
                                ##portion of the names in the bands
                                ##must match each other. Names of the
                                ##bands should also contain the band
                                ##numbers. Only bands with similar
                                ##extents are stacked while ignoring
                                ##data with different
                                ##extents. \code{\link{equateEEnames}}
                                ##is implemented to match data from
                                ##different satellites
(
    tarPath, ##<<\code{List}. File paths or names of the compressed
             ##(\code{.tar}) data containing multi- or hyper-spectral
             ##radiometric data.
    bands = NULL, ##<<\code{numeric} or \code{NULL}. Bands to be
                  ##processed after matching the band names (see
                  ##\code{sat.nm}). Default \code{NULL} preserves the
                  ##bands.
    sat.nm = 'LE07',  ##<<\code{character} or \code{NULL}. Regular
                      ##expression indicating a Landsat satellite used
                      ##to match names of the bands. Default
                      ##\code{'LE07'}. If \code{NULL} then the band
                      ##names are not matched.
    drop.unmatched = TRUE, ##<<\code{logical}. Drop the unmatched
                          ##bands.
    meta.wd = FALSE   ##<<\code{logical}. Copy metadata in the working
                     ##directory.

){
    fimp <- function(tarPath, bands){
        if(is.null(tarPath))
            return(NULL)
        tarNm <- basename(tarPath)
        ## tmp. <- tempdir()
        ## temp. <- file.path(tmp., 'EEstack')
        lstar <- untar(tarPath, exdir = temp.,
                       list = TRUE)
        if(!all(lstar%in%list.files(temp.)))
            untar(tarPath, exdir = temp.)
        cmmPth <- Reduce(PTXQC::LCS, lstar)
        tif <- lstar[grepl('.tif', lstar,
                           ignore.case = TRUE)]

        ftiford <- function(tif){
            tifa <- tif
            cmmPth <- Reduce(PTXQC::LCS, tif)
            tif <- gsub(cmmPth,'', tif)
            names(tif) <- 1:length(tif)
            fOrd <- function(x)
                x[order(nchar(x), x)]
            tif. <- fOrd(tif[grepl('[0-9]', tif)])
            tif.. <- fOrd(tif[!tif%in%tif.])
            tif.. <- as.numeric(names(c(tif., tif..)))
            return(tifa[tif..])}
        tif <- ftiford(tif)

        pathtifs <- file.path(temp., tif)
        rasterImp <- Map(function(x)
            raster::raster(x), pathtifs)
        exts <- Map(extent, rasterImp)
        isext <- sapply(exts, function(x)
            extent(rasterImp[[1L]]) == x)
        rasterImp <- rasterImp[isext]
        tifNm. <- gsub(cmmPth,'', tif)
        tifNm <- gsub('^_|.tif','', tifNm.,
                      ignore.case = TRUE)
        names(rasterImp) <- tifNm[isext]
        stkk <- stack(rasterImp)
        attributes(stkk) <- c(attributes(stkk), list(nam.= cmmPth))
        ## unlink(temp., recursive = TRUE)
        return(stkk)
    }
    tmp. <- tempdir()
    temp. <- file.path(tmp., 'EEstack')

    mta <- Map(function(x,y)
        fimp(x,y), tarPath,
        MoreArgs = list(y = bands))
    mtt <- lapply(mta, function(x)
        attr(x, 'nam.'))
    names(mta) <- unlist(mtt)
    nm. <- names(mta)
    if(!is.null(sat.nm))
        mta <- equateEEnames(mta, sat.nm)

    excNas <- function(x){
        nms. <- names(x)
        nam. <- !grepl('NA', nms.)
        nbds <- raster::subset(x, nms.[nam.])
        return(nbds)}
    if(drop.unmatched)
    mta <- Map(function(x)
        excNas(x),as.list(mta))
    if(meta.wd){
        lstf <- list.files(temp.,'MTL.txt')
        fils <- file.path(temp., lstf)
        file.copy(fils, getwd())
    }
    if(!is.null(bands))
        mta <- Map(function(x)x[[bands]],
                   mta)
    class(mta) <- c("EEstack",class(mta))
    return(mta)
### RasterBrick.
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
})
