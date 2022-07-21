


#' Sample from uniform distribution by optim
#'
#' @param xvec A vector of min,max values . Ex:c(10,20)
#' @param N Desired number of values
#' @param medValue  Desired Median 
#'
#' @return Samples with approx desired median. 

optimUnisampling <- function(xvec,N=1000,medValue,floor=TRUE)
{
  
  #Function to calculate distance between sampling median and desired
  fun <- function(xvec, N=1000) {
    xmin <- xvec[1]
    xmax <- xvec[2]
    if (floor){
    x <- floor(runif(N, xmin, xmax))}
    else{
      x <- runif(N, xmin, xmax) 
    }
    xdist <- (median(x)-medValue)^2
    xdist
  }
# Optimization 
xr <- optim(xvec, fun) 
xrmin <- xr$par[[1]]
xrmax <- xr$par[[2]]
sampled <- runif(N, min = xr$par[[1]], max = xr$par[[2]])
if (xrmin==xvec[1] & xrmax==xvec[2] & floor)    return (floor(sampled))
else if (xrmin==xvec[1] & xrmax==xvec[2])   return (sampled)
else return (optimUnisampling(xvec,N=1000,medValue))
}


### Sample Weight 

#wt <- optimUnisampling(xvec=c(31,157),N=1000,medValue = 75)
#age <- optimUnisampling(xvec=c(19,90),N=1000,medValue = 74)
#height <- optimUnisampling(xvec=c(122,201),N=1000,medValue = 168)
#bmi <- optimUnisampling(xvec=c(13.4,56.9),N=1000,medValue = 25.9,floor=FALSE)
#creatCl <- optimUnisampling(xvec=c(6.7,140),N=1000,medValue = 71.7,floor=FALSE)
#sBil <- optimUnisampling(xvec=c(0.1,2.4),N=1000,medValue = 0.5,floor=FALSE)
#sAlb <- optimUnisampling(xvec=c(1.9,5.0),N=1000,medValue = 4.0,floor=FALSE)
#sAlat <- optimUnisampling(xvec=c(5,287),N=1000,medValue = 19,floor=FALSE)




