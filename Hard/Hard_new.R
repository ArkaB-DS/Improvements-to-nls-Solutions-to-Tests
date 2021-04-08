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
lwd=2,col="red",type="o",xlim=c(1,37))

# clean the data
time <- time[y>0]
y <- y[y>0]

## model based on a list of parameters
getPred <- function(parS, x) {
 return ( parS$a /(1+  parS$b*exp(-x * parS$c)) )
}

## residual function
residFun <- function(p, observed, x) {
observed- getPred(p,x)
}

## starting values for parameters
parStart <- list(a=10.04212,b=1352.3,c=0.33883)

## fit using internal approximation(no Jacobian)
model<- y ~ a/(1+b*exp(-c*time))
Data=data.frame(y,time)
nlmod <- try(nlxb(model,
start=parStart, trace=FALSE,data=Data))
print(nlmod)
y_pred<-predict(nlmod,Data)
lines(y_pred,col="blue",pch=16,lwd=2)

# to find the Jacobian function
Z<-model2rjfun(modelformula=model, pvec=parStart, data = Data, jacobian = TRUE,
 testresult = TRUE)
jfn<-modelexpr(Z)

# the jacobian function is- 
jacobian<-function(x,observed,Pars){
    .expr3 <- exp(-Pars$c * x)
    .expr5 <- 1 + Pars$b * .expr3
    .expr10 <- .expr5^2
    .value <- Pars$a/.expr5 - observed
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a", 
        "b", "c")))
    .grad[, "a"] <- 1/.expr5
    .grad[, "b"] <- -(Pars$a * .expr3/.expr10)
    .grad[, "c"] <- Pars$a * (Pars$b * (.expr3 * x))/.expr10
    .grad
}

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
