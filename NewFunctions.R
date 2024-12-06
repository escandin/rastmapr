int2bits=function (x, n, nbits=16)
  # Converts an integer into the corresponding bit denomination for a number of bits specified by n
  # x: integer to be converted into bits
  # n: bit numbers to be retreived. It can be a number or a vector specifying the bit positions to retrieve.
  #   The first bit is counted as zero and they are read from right to left
{
  bit <-intToBits(x)
  rev(as.numeric(paste(tail(rev(as.integer(bit)), nbits))))[n+1]
  #if (reverse==TRUE){bit=rev(bit)}
  #return(bit)
}

bitMask=function(r, bit){
  # produces a mask that excludes as NA, pixels with a value of 1 in the bit position 
  # designated in the argument "bit" in the landsat QA layer.
  # All other pixels are assigned a value of 1.
  x <- calc(r, fun = function(x) int2bits(x, bit))
  v <- raster::getValues(x)
  v[v==1]=NA
  v[v==0]=1
  x <- raster::setValues(x, v)
  return(x)
}

bitMask2=function(r, bit){
  # produces a mask that excludes as NA, pixels with a value of 1 in the bit position 
  # designated in the argument "bit" in the landsat QA layer
  x <- raster::getValues(x)
  v=data.frame(lapply(x,function(x) int2bits(x, bit)))
  v[v==1]=NA
  v[v==0]=1
  x <- raster::setValues(x, v)
  return(x)
}

test=raster::getValues(L8rsz[[8]])
out=data.frame(lapply(test,function(x) int2bits(x, 4)))
test=setValues(as.vector(test), out)


