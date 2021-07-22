equateEEnames <- structure(function #Equate EE names
### This function equates names of bands in Earth-Explorer (EE) data
### sets using conventions of a specific Landsat satellite.
                           ##details<< Landsat satellite in
                           ##\code{sat.nm} is matched with columns in
                           ##\code{\link{LBands}}. 
(
    dtls, ##<<\code{list}. Set of two \code{RasterBrick} objects:
            ##target and reference, sharing a common region such as
            ##these produced by \code{\link{EEstackWithoutMeta}} or
            ##\code{\link{RasterIntersection}}.
    sat.nm = 'LE07', ##<<\code{character}. Regular expression
                    ##indicating a Landsat satellite used as reference
                    ##to match names of the bands, see
                    ##\code{Details}.
    change.dt = TRUE ##<<\code{logical}. Print the renamed data. If
                     ##\code{FALSE} then only new names are printed.
){
    ## mtptt <- substr(names(dtls),1,4)
    ## if(Reduce('==', mtptt)){
    ##     sat.nm <- mtptt[1L]
    ##     print('equal satellite; preserving names')
    ## }
    equate <- function(dtls, pos. = 1L, sat.nm, change.dt){
    data('LBands', envir = environment())
    dtls <- dtls[pos.]
    sat.nm. <- grepl(sat.nm, get('LBands')$'patt')
    dbref <- subset(get('LBands'), sat.nm.)
    dbother <- subset(get('LBands'), !sat.nm.)
    dbmr <- merge(dbother, dbref[, c('name','band')],
                  by = 'name', all = TRUE)
    dbmr <- dbmr[order(dbmr$'patt', dbmr$'band.x'),]
    mtptt <- substr(names(dtls),1,4)
    if(mtptt == sat.nm)
        if(change.dt){
            return(dtls)
        } else {
            return(names(dtls[[1L]]))}
    cases <- mapply(function(x)
        grepl(x, mtptt), dbmr$'patt',
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE)
    dfr <- subset(dbmr[, c('band.x','band.y')], cases%in%1)
    nm.. <- names(dtls[[1L]])
    np = '[0-9]'
    nmn <- nm..[grepl(np, nm..)]
    nmc <- nm..[!grepl(np, nm..)]
    extn <- as.numeric(gsub("[^0-9]", "", nmn))
    subsn <- gsub(np, "", nmn)
    neon. <- dfr[dfr$'band.x'%in%extn, 'band.y']
    newalph <- paste(subsn, neon., sep = '')
    newnam <- c(newalph, nmc)
    if(change.dt){
        names(dtls[[1L]]) <- newnam
        return(dtls)}
    return(newnam)}
    
    eqt <- list()
    for(i in 1: length(dtls))
        eqt[i] <- equate(dtls, pos. = i, sat.nm, change.dt)
    names(eqt) <- names(dtls)
    return(eqt)
### \code{list} of either characters or \code{RasterBrick} depending
### on \code{change.dt}
} , ex=function(){
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
                  'LT050060661988072201T1-SC20181031160603.tar.gz')
    tarPaths <- system.file(tarFiles, package = 'aRn')
    stack <- EEstackWithoutMeta(tarPaths, c(1,3:6), sat.nm = NULL)
    eqstack <- equateEEnames(stack)
})
