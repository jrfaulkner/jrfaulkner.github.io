# LatticeKrig examples
library(LatticeKrig)
library(fields)
library(latex2exp)

## example 1
# Load ozone data set
data(ozone2)  
x<-ozone2$lon.lat
y<- ozone2$y[16,]
# Find location that are not 'NA'.
# (LKrig is not set up to handle missing observations.)
good <-  !is.na( y)
x<- x[good,]
y<- y[good]

quilt.plot(x, y)
US(add=TRUE)

# thin plate spline-like model with the lambda parameter estimated by
# maximum likelihood. Default choices are made for a.wght, nlevel, NC
# and alpha.

obj<- LatticeKrig( x, y)
## Not run: 
# summary of fit and a plot of fitted surface
print( obj)
surface( obj )
US(add=TRUE)
points(x)

# prediction standard errors at observations
out.se<- predictSE( obj, xnew= x)
# predict mean at observations:
out.fhat<- predict( obj, xnew= x)
# coveniently predict on a 100X100 grid for plotting
out.surf<- predictSurface( obj, nx=100, ny=100)
surface(out.surf)

## example 2:
LKinfo<- LKrigSetup( x, nlevel=3, nu=1, NC=5, a.wght=5,
                     lambda=.01)
# maximize likelihood over lambda see help( LKrig.MLE) for details
# this assumes the value of 5 for a.wght.  In many cases the fit is not
# very sensitive to the range parameter such as a.wght in this case --
# but very sensitive to lambda when varied on a log scale.

MLE.fit<- LKrig.MLE(x,y, LKinfo=LKinfo)
MLE.fit$summary # summary of optimization over lambda.
# fit using MLE for lambda MLE function has added MLE value of lambda to
# the LKinfo object.
obj<- LKrig( x,y, LKinfo=MLE.fit$LKinfo.MLE)  
print( obj) 
out.surf<- predictSurface( obj, nx=100, ny=100)
surface(out.surf)
US(add=TRUE)
points(x)

## example 3:

########################################################################
# A bigger problem: 26K observations and 4.6K basis functions
# fitting takes about 15 seconds on a laptop for a fixed covariance
#  LKrig.MLE to find the MLE (not included) for lambda takes abou
#  8 minutes
#######################################################################
## Not run: 
data(CO2)
# the Ring geometry will be periodic in the first dimension and rectagular on 
# second. A useful approximation for spherical data omitting the polar caps. 

LKinfo.CO2<- LKrigSetup(CO2$lon.lat, NC=100,nlevel=1, lambda=.2,
                        a.wght=5, alpha=1, 
                        LKGeometry="LKRing", choleskyMemory = list(nnzR=2e6) )
print(LKinfo.CO2)                                          
obj1<- LKrig( CO2$lon.lat,CO2$y,LKinfo=LKinfo.CO2)
# 5700+ basis functions 101X57 lattice  (number of basis functions
# reduced in y direction because of a rectangular domain
obj1$trA.est # about 2900+ effective degrees of freedom 
#
glist<- list( x= seq( -180,180,,200),y=seq( -80,80,,100) )
fhat<- predictSurface( obj1,grid.list=glist)
#Plot data and gap-filled estimate
set.panel(2,1)
quilt.plot(CO2$lon.lat,CO2$y,zlim=c(373,381))
title("Simulated CO2 satellite observations")
world(add=TRUE,col="magenta")
image.plot( fhat,zlim=c(373,381))
world( add=TRUE, col="magenta")
title("Gap-filled global predictions")



### other plots
set.panel(1,1)
xs = seq(0, 1.3, l=100)
wend = function(x) {
  out = (1-x)^6 * (35*x^2 + 18*x + 3) / 3
  inDom = (x >= 0) & (x <= 1)
  out[!inDom] = 0
  out
}
plot(xs, wend(xs), main="Wendland Function", xlab="d", ylab=TeX("$\\phi(d)$"), type="l")
abline(v=1, col="red", lty=2)
