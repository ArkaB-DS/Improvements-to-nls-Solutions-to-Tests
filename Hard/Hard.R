### "Improvements to nls()" - GSOC'21
## Hard

# Install and load necessary packages
install.packages("nlsr")
install.packages("minpack.lm")
library(nlsr)
library(minpack.lm)

# test data

time <- c( 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
 30, 31, 32, 33, 34, 35)
y <- c( 0.0074203,  0.3188325,  0.2815891, -0.3171173, -0.0305409,  0.2266773,
  -0.0216102,  0.2319695, -0.1082007,  0.2246899,  0.6144181,  1.1655192,
  1.8038330,  2.7644418,  4.1104270,  5.0470456,  6.1896092,  6.4128618,
  7.2974793,  7.8965245,  8.4364991,  8.8252770,  8.9836204, 9.6607736,
  9.1746182,  9.5348823, 10.0421165,  9.8477874,  9.2886090,  9.3169916,
  9.6270209 )


# Plot y vs. time
#par(mfrow=c(1,2))
plot(time,y,xlab="t",ylab=expression(y[t]),
main="Plot of observed and fitted y vs. Time",
lwd=2,col="red",type="o",xlim=c(1,37))

# clean the data
y.max=max(y) # 10.04212
time <- time[y>=0&y!=max(y)]
y <- y[y>=0&y!=max(y)]


lmod<-lm( log( (y/y.max) / (1-(y/y.max)) ) ~ time ) 

## model based on a list of parameters
getPred <- function(parS, x) {
 return ( parS$a /(1+  parS$b*exp(-x * parS$c)) )
}

## residual function
residFun <- function(p, observed, x) {
observed - getPred(p,x)
}

## jacobian function
jbn<-function(par,x,observed){
return( c(-1/(1+par$b*exp(-x*par$c)),
	par$a*exp(-x*par$c)/(1+par$b*exp(-x*par$c))^2,
	-par$a*par$b*x*exp(-x*par$c)/(1+par$b*exp(-x*par$c))^2
) )
}


## starting values for parameters
parStart <- list(a=10.04212,b=exp(6.7599), c=0.3147)

nls.out <- nls.lm(par=parStart, fn = residFun, observed = y,
x = time, control = nls.lm.control(nprint=1),jac=jbn)
summary(nls.out)
nls.out$par

# get the fitted values of y
fit.y<-getPred(nls.out$par,time)
lines(fit.y,col="blue",pch=16,lwd=2)
legend("topleft",legend=c("observed y","fitted y"),
pch=c(16,NA),lwd=2,col=c("red","blue"))

#modelexpr(model2rjfun(log( (y/y.max) / (1-(y/y.max)) ) ~ time,c(a,b,c),
#jacobian=T))

#nls.out <- nlfb(parStart, resfn = residFun, data=data.frame(y,time),
# control = nls.lm.control(nprint=1),jacfn=jbn)
#summary(nls.out)
nls.out$par