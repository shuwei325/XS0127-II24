
# Laboratorio 5 -----------------------------------------------------------


library(ggplot2)
library(forecast)
library(fpp2)

# 1. Ejemplo de ruido blanco ---------------------------------------------

w = rnorm(500,0,1) 
ts.plot(w)
acf(w,lag.max = 50,main="estimación de acf de w")


# 2. Ejemplo de medias móviles --------------------------------------------

w = rnorm(500,0,1) # 500 N(0,1) variates
v = stats::filter(w, sides=2, filter=rep(1/3,3)) # moving average
v = na.omit(v)
ts.plot(v)
acf(v,lag.max = 50,main="estimación de acf de v")

w = rnorm(500,0,1) # 500 N(0,1) variates
v = filter(w, sides=2, filter=rep(1/7,7)) # moving average
v = na.omit(v)
ts.plot(v)
acf(v,lag.max = 50,main="estimación de acf de v")



# 8. Ruido blanco y las medias móviles -------------------------------------------------------------

w = rnorm(500,0,1) # 500 N(0,1) variates
v = stats::filter(w, sides=2, filter=rep(1/3,3)) # moving average
par(mfrow=c(2,1))
plot.ts(w, main="white noise")
plot.ts(v, ylim=c(-3,3), main="moving average")

#las funciones de autocorrelación
acf(w)
acf(na.omit(v))

plot.ts(AP.ts)
AP.v1 = stats::filter(AP.ts, sides=2, filter=rep(1/3,3)) # moving average
AP.v2 = stats::filter(AP.ts, sides=2, filter=rep(1/6,6)) # moving average
AP.v3 = stats::filter(AP.ts, sides=2, filter=rep(1/12,12)) # moving average
points(AP.v1,type="l",col=2)
points(AP.v2,type="l",col=3)
points(AP.v3,type="l",col=4)
legend("topleft",legend=c("MA-3","MA-6","MA-12"),
       col=c(2,3,4),lty=1)



# 9. Señal+ruido ----------------------------------------------------------

cs = 2*cos(2*pi*1:500/50 + .6*pi); w = rnorm(500,0,1)
par(mfrow=c(3,1), mar=c(3,2,2,1), cex.main=1.5)
plot.ts(cs, main=expression(2*cos(2*pi*t/50+.6*pi)))
plot.ts(cs+w, main=expression(2*cos(2*pi*t/50+.6*pi) + N(0,1)))
plot.ts(cs+5*w, main=expression(2*cos(2*pi*t/50+.6*pi) + N(0,25)))




# 3. Ejemplo de sorteo navideño -------------------------------------------

sorteo<-read.csv("sorteo.csv",sep=",")
y<-ts(sorteo$numero)
autoplot(y)
ggAcf(y,lag.max = 50,main="estimación de acf")


# 4. Ejemplo de graduados de ITCR -----------------------------------------

itcrgrad<-read.csv("ITCR.csv",sep=",")
y<-ts(itcrgrad$graduados,start=1975)
autoplot(y) 
ggAcf(y)
w=diff(y)
autoplot(w) 
ggAcf(w)

# 5. Ejemplo de turistas -----------------------------------------

turistas<-read.csv("turistas.csv",sep=";")
y<-ts(turistas$turistas,start=c(1991,1),frequency=12)
autoplot(y) 
ggAcf(y)
w=diff(y)
autoplot(w) 
ggAcf(w)

autoplot(log(y)) 
ggAcf(log(y))
logw=diff(log(y))
autoplot(logw) 
ggAcf(logw,lag.max=50)


