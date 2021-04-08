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
plot(time,y,xlab="t",ylab=expression(y[t]),
main="Plot of observed and fitted y vs. Time",
lwd=2,col="red",type="o",xlim=c(5,35))

# clean the data
time <- time[y>0]
y <- y[y>0]

## model based on a list of parameters
getPred <- function(parS, x) {
  return ( parS$a /(1+  parS$b*exp(-x * parS$c)) )
}

## residual function
residFun <- function(p, observed, x) {
  return(observed- getPred(p,x))
}

## starting values for parameters
parStart <- list(a=10.04212,b=1352.3,c=0.33883)

## 
model<- y ~ a/(1+b*exp(-c*time))
Data=data.frame(time,y)

# the analytic jacobian function is- 
jacobian<-function(x,observed,Pars){
	mat <- matrix(0,nrow=length(x),ncol=length(Pars))
      colnames(mat)<-c("a","b","c")
	mat[,"a"]<- -1/(1+Pars$b*exp(-Pars$c*x))
	mat[,"b"]<- Pars$a*exp(-Pars$c*x)/(1+Pars$b*exp(-Pars$c*x))^2
	mat[,"c"]<- -Pars$a*x*Pars$b*exp(-Pars$c*x)/(1+Pars$b*exp(-Pars$c*x))^2
      return(mat)
}

# the approximate jacobian function is- 
jacobian.approx<-function(x,observed,parS){
	delta<- 5
	mat <- matrix(0,nrow=length(x),ncol=length(parS))
      colnames(mat)<-c("a","b","c")
	mat[,"a"]<-  ((parS$a+delta) /(1+  parS$b*exp(-x * parS$c))-
			 parS$a /(1+  parS$b*exp(-x * parS$c)))/delta
	mat[,"b"]<- (parS$a /(1+  (parS$b+delta)*exp(-x * parS$c))-
			 parS$a /(1+  parS$b*exp(-x * parS$c)))/delta
	mat[,"c"]<- (parS$a /(1+  parS$b*exp(-x * (parS$c+delta)))-
			 parS$a /(1+  parS$b*exp(-x * parS$c)))/delta
      return(mat)
}

# fitting the nonlinear model based on approximate jacobian
nls.out.approx <- nls.lm(par=parStart, fn = residFun, observed = y,
x = time, control = nls.lm.control(nprint=1),jac=jacobian.approx)
summary(nls.out.approx)
nls.out.approx$par
a.est.approx=nls.out.approx$par[["a"]]
b.est.approx=nls.out.approx$par[["b"]]
c.est.approx=nls.out.approx$par[["c"]]
lines(time,a.est.approx/(1 + b.est.approx * exp(-c.est.approx * time)),
col="blue",pch=16,lwd=2)

# fitting the nonlinear model based on Jacobian
nls.out <- nls.lm(par=parStart, fn = residFun, observed = y,
x = time, control = nls.lm.control(nprint=1),jac=jacobian)
summary(nls.out)
nls.out$par
a.est=nls.out$par[["a"]]
b.est=nls.out$par[["b"]]
c.est=nls.out$par[["c"]]
lines(time,a.est/(1 + b.est * exp(-c.est * time)),col="seagreen",pch=16,lwd=2)
legend("topleft",legend=c("observed y","fitted y(approx)","fitted y(Jacobian)"),
pch=c(16,NA,NA),lwd=2,col=c("red","blue","seagreen"),cex=0.8)
