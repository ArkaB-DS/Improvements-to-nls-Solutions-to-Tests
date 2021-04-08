### "Improvements to nls()" - GSOC'21
## Medium

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

# Plot observed y vs. time
plot(time,y,xlab="t",ylab=expression(y[t]),
main="Plot of Observed and Fitted y against Time",
lwd=2,col="red",type="o")

# clean the data
time <- time[y>=0]
y <- y[y>=0]

# Approximate estimation of a,b,c 

# approximate value of a
a<-max(y)

# approximate value of b
b <- a/y[1]-1

# approximate value of c

time<-time[y!=a]
y<-y[y!=a]

c<-log( (b*y)/(a-y) )/time
#c<-mean(c[-length(c)])
c<-mean(c)
# plot of fitted y values
lines(time, a/(1+b*exp(-c*time)),type="o",lwd=2,col="blue")
legend("topleft",legend=c("Observed y","Fitted y"),
col=c("red","blue"),lwd=2,,pch=16)

# print a,b,c
print(a)
print(b)
print(c)
