rclMatrix <- structure(function #reclassification matrix
### A binary reclassification matrix is produced by specifying a
### constant (threshold). The constant is used to split the interval
### \code{c(-Inf, Inf)} into two intervals.
                                ##details<< See \code{rcl} argument in
                                ##\code{\link{reclassify}}. This
                                ##function is implemented by
                                ##in-package \code{\link{PIFmodel}}
                                ##and \code{\link{RasterIntersection}}
(
    thr = 0, ##<<\code{numeric}. Constant (threshold). Default \code{0}
    oneFirst = TRUE ##<<\code{numeric}. Condition in the first row
              ##becomes \code{1}. If \code{FALSE} then the second
              ##column will do.
){
    m <- diag(c(-1,1) * Inf)
    then <- c(1, NA)
    if(thr != 0)
        m[m == 0] <- thr
    if(!oneFirst)
        then <- rev(then)
    mat <- cbind(m, then)
    nam. <- c('from','to','becomes')
    colnames(mat) <- nam.
    return(mat)
### \code{Matrix}.
} , ex=function(){
    rclMatrix(4, FALSE)
})
