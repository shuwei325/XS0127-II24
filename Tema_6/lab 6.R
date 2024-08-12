
# Laboratorio 6 -----------------------------------------------------------

library(ggplot2)
library(forecast)
library(fpp2)
library(astsa)
library(car)
library(TSA)
library(tseries)
library(urca)

# 1. Modelo ARIMA -------------------------------------------------

# 1.1 ejemplo simulado ----------------------------------------------------

# generación de datos
gen_ar1a <- function(N = 150, phi1 = 0.8, sigma2 = 1) {
  a <- rnorm(N,0,sigma2) 
  y <- as.numeric(0)
  y[1] <- a[1]
  for(i in 2:N){
    y[i] <- phi1*y[i-1]+a[i]
  }
  return(y)
}                     

gen_ar1b <- function(N = 150, C=0, phi1 = 0.8, sigma2 = 1) {
  NN <- 1000
  a <- rnorm(NN+N,0,sigma2) 
  y <- as.numeric(0)

  y[1] <- a[1]
  for(i in 2:(NN+N)){
    y[i] <- C + phi1*y[i-1] + a[i]
  }
  return(y[NN:(NN+N)])
}

phi1=0.6

y <- gen_ar1b(N=150,C=5,phi1=phi1,sigma2=1)

# descriptiva
ts.plot(y)
mean(y)  #promedio teórico: 5/(1-phi1)
acf(y,lag.max=30)
pacf(y)
acf2(y)   #library(astsa)

?Arima #library(forecast)
mod0a <- Arima(y, order=c(1,0,0),method="CSS-ML")
mod0b <- Arima(y, order=c(1,0,0),method="ML")
mod0c <- Arima(y, order=c(1,0,0),method="CSS")
summary(mod0a)
summary(mod0b)
summary(mod0c)

#Note que no está estimando C!
mod0d <- Arima(y, order=c(1,0,0))
summary(mod0d)

mean(y)
5/(1-phi1) #media teórica

?tseries::arma
mod0e<-tseries::arma(y,order=c(1,0),include.intercept=TRUE)
summary(mod0e)

#devolvemos al mod0d
res<-mod0d$residuals
ts.plot(res)
acf(res)
pacf(res)
acf2(res)

tsdiag(mod0d)    #library(stats)
checkresiduals(mod0d,lag=10)
checkresiduals(mod0d,lag=30)
#Contraste de Ljung-Box

#Normalidad
shapiro.test(res)
jarque.bera.test(res)
qqPlot(res)

#Pronóstico
forecast(mod0d)
autoplot(forecast(mod0d))


# 1.2 ejemplo simulado con R ---------------------------------------------------

#AR(1)
m<-5 #la media del proceso
y1 <- arima.sim(n = 150, model = list(order = c(1,0,0),ar = c(0.8)),sd=3,rand.gen= rnorm) + m
ts.plot(y1)
acf2(y1)

mod1<- forecast::Arima(y1, order = c(1, 0, 0))
summary(mod1)


#AR(2) #ejercicio tomado de Bernhard (2008)
series <-  rnorm(1000)
y.st <- filter(series, filter=c(0.6, -0.28),
               method='recursive')
ar2.st <- arima(y.st, c(2, 0, 0), include.mean=FALSE,
                transform.pars=FALSE, method="ML")
ar2.st$coef
polyroot(c(1, -ar2.st$coef))
Mod(polyroot(c(1, -ar2.st$coef)))
root.comp <- Im(polyroot(c(1, -ar2.st$coef)))
root.real <- Re(polyroot(c(1, -ar2.st$coef)))
# Plotting the roots in a unit circle
x <- seq(-1, 1, length = 1000)
y1 <- sqrt(1- x^2)
y2 <- -sqrt(1- x^2)
plot(c(x, x), c(y1, y2), xlab='Real part',
     ylab='Complex part', type='l',
     main='Unit Circle', ylim=c(-2, 2), xlim=c(-2, 2))
abline(h=0)
abline(v=0)
points(Re(polyroot(c(1, -ar2.st$coef))),
       Im(polyroot(c(1, -ar2.st$coef))), pch=19)
legend(-1.5, -1.5, legend="Roots of AR(2)", pch=19)

#Otra posibilidad es usar el inverso de las raíces.
autoplot(ar2.st)


#ARMA(1,1)
y2<-arima.sim(n = 150, list(order = c(1,0,1),ar = c(0.88), ma = c(-0.23)),
              sd = sqrt(2))

ts.plot(y2)
acf(y2)
pacf(y2)
acf2(y2)

#ARMA(1,1)
mod2a<- forecast::Arima(y2, order = c(1, 0, 1))
summary(mod2a)
checkresiduals(mod2a,lag=10)

#AR(1)
mod2b<- forecast::Arima(y2, order = c(1, 0, 0))
summary(mod2b)
checkresiduals(mod2b,lag=10)

mod2a$aic
mod2b$aic

#procedimiento automático (pero tener mucho cuidado!!!)
auto.arima(y2,ic="aicc") #por defecto
auto.arima(y2,ic="aic")
auto.arima(y2,ic="bic")

###Contraste de raíz unitaria
#Probamos con dos tamaño de series
TT=150
TT=500

#AR(1)
y1 <- arima.sim(n = TT, model = list(order = c(1,0,0),ar = c(0.8)),sd=3,rand.gen= rnorm)
ts.plot(y1)
adf.test(y1) 

#ARIMA(0,1,0)
y2 <- arima.sim(n = TT, model = list(order = c(0,1,0),sd=1,rand.gen= rnorm))
ts.plot(y2)
acf2(y2)
adf.test(y2)


# 2. Ejemplo con graduados de ITCR de 1975 a 2002 -----------------------


itcrgrad<-read.csv("ITCR.csv",sep=",") 
y<-ts(itcrgrad$graduados,start=1975)

ts.plot(y) 
acf(y)
pacf(y)
acf2(y)
#indicación de no estacionariedad. 
# Como ejemplo vamos a ajustar un AR(1) (1 rezago de f.a.c.p. significativo)

mod0 <- Arima(y, order=c(1,0,0))

mod0a <- Arima(y, order=c(1,0,0),method="CSS-ML")
mod0b <- Arima(y, order=c(1,0,0),method="ML")
mod0c <- Arima(y, order=c(1,0,0),method="CSS")
summary(mod0b)
summary(mod0c)

adf.test(y)

dif.y<-diff(y)

ts.plot(dif.y)
acf(dif.y)
pacf(dif.y)

adf.test(dif.y)

mod1 <- Arima(y, order=c(0,1,0))
summary(mod1)

res<-mod1$res
ts.plot(res)
acf(res)
pacf(res)


tsdiag(mod1)
checkresiduals(mod1,lag=10)

#Normalidad
shapiro.test(res)
jarque.bera.test(res)
qqPlot(res)

#Pronóstico
forecast(mod1)
autoplot(forecast(mod1))


#procedimiento automático
auto.arima(y,ic="aicc") #por defecto
auto.arima(y,ic="aic")
auto.arima(y,ic="bic")

auto.arima(y,ic="aicc", allowdrift = FALSE) #por defecto
auto.arima(y,ic="aic", allowdrift = FALSE)
auto.arima(y,ic="bic", allowdrift = FALSE)

# 3. Tasa de desempleo --------------------------------------------------

# Ejemplo 1-3 tomado de Bernhard (2008): Tasa de desempleo

data(npext)
y <- ts(na.omit(npext$unemploy), start=1909, end=1988,
        frequency=1)
plot(y, ylab="unemployment rate (logarithm)")
acf2(y,ylim=c(-1, 1))

#1.  ARMA(2,0) 
arma20 <- Arima(y, order=c(2, 0, 0))
summary(arma20)

loglik <- arma20$loglik
(aic<- -2*loglik+2*(2+1+1))
arma20$aic

TT <- length(y)
(aicc <- aic+(2*(2+1+1)*(2+1+2))/(TT-2-1-2))
arma20$aicc

(bic <- aic+(log(TT)-2)*(2+1+1))
arma20$bic

res20 <- residuals(arma20)
ts.plot(res20)
shapiro.test(res20)
qqPlot(res)
tsdiag(arma20)
checkresiduals(arma20)
arma20$coef
autoplot(arma20)


#2.  ARMA(1,1)
arma11 <- Arima(y, order = c(1, 0, 1))
summary(arma11)
tsdiag(arma11)

c(arma20$aic,arma20$aicc,arma20$bic)
c(arma11$aic,arma11$aicc,arma11$bic)


res11 <- residuals(arma11)
ts.plot(res11)
shapiro.test(res11)
tsdiag(arma11)
checkresiduals(arma11)
autoplot(arma11)

arma20$aic
arma11$aic

## auto.arima()
arma.auto<-auto.arima(y, max.p = 3, max.q = 3, start.p = 1,
                      start.q = 1, ic = "aic")

arma.auto

## Forecasts
arma11.pred <- predict(arma11, n.ahead = 10)
predict <- ts(c(rep(NA, length(y) - 1), y[length(y)],
                arma11.pred$pred), start = 1909,
              frequency = 1)
upper <- ts(c(rep(NA, length(y) - 1), y[length(y)],
              arma11.pred$pred + 2 * arma11.pred$se),
            start = 1909, frequency = 1)
lower <- ts(c(rep(NA, length(y) - 1), y[length(y)],
              arma11.pred$pred - 2 * arma11.pred$se),
            start = 1909, frequency = 1)
observed <- ts(c(y, rep(NA, 10)), start=1909,
               frequency = 1)
## Plot of actual and forecasted values
plot(observed, type = "l",
     ylab = "Actual and predicted values", xlab = "")
lines(predict, col = "blue", lty = 2)
lines(lower, col = "red", lty = 5)
lines(upper, col = "red", lty = 5)
abline(v = 1988, col = "gray", lty = 3)

plot(forecast(arma11))



# 4. Producto nacional bruto - U.S. (ejemplo 3.40, Shumway&Stoffer) ----------------------------------------------------------------
# producto nacional bruto, U.S. (en mil millones y son datos trimestrales de 1947 a 2002)
# los datos fueron ajustada estacionalmente.

y<-astsa::gnp

ts.plot(y)
acf2(y, 50)           

#contraste de Dickey-Fuller
adf.test(y)
ts.plot(diff(y))
ts.plot(log(y))
adf.test(log(y))


dif.log.y = diff(log(y))      # growth rate
plot(dif.log.y)
acf2(dif.log.y, 24)  

#Contraste de Dickey-Fuller
adf.test(dif.log.y)

#AR(1)
moda<-Arima(dif.log.y, order=c(1,0,0))
summary(moda)
autoplot(moda)

#MA(2)
modb<-Arima(dif.log.y, order=c(0,0,2))
summary(modb)
autoplot(modb)

checkresiduals(moda)
checkresiduals(modb)

c(moda$aic,moda$aicc,moda$bic)
c(modb$aic,modb$aicc,modb$bic)


modc<-Arima(log(y),order=c(1,1,0),include.drift=TRUE)
summary(modc)

modd<-Arima(y,order=c(1,1,0),include.drift=TRUE,lambda=0)
summary(modd)

mode<-Arima(y,order=c(0,1,2),include.drift=TRUE,lambda=0)
summary(mode)


#Modd
checkresiduals(modd,lag=10)
checkresiduals(mode,lag=10)

res<-modd$res
#res<-mode$res

ts.plot(res)
acf(res)
pacf(res)

#Normalidad
shapiro.test(res)
jarque.bera.test(res)
qqPlot(res)

#¿Cuál modelo es mejor?


# el parquete astsa
astsa::sarima(dif.log.y, 1, 0, 0)
astsa::sarima(dif.log.y, 0,0, 2)

#Pronóstico
forecast(moda)
autoplot(forecast(moda))

forecast(modd)
autoplot(forecast(modd))


# 5.logarítmo de varve glacial   (ejemplo 3.41, Shumway&Stoffer)-----------------------------------------------
# Los glaciares que se derriten depositan capas anuales de arena y limo durante 
# las temporadas de derretimiento de primavera, que pueden reconstruirse anualmente 
# durante un período que va desde el momento en que comenzó la desglaciación en Nueva 
# Inglaterra (hace unos 12.600 años) hasta el momento en que terminó (hace unos 6000 años).
# Dichos depósitos sedimentarios, llamados varvas, pueden utilizarse como sustitutos de 
# parámetros paleoclimáticos, como la temperatura, porque, en un año cálido, se depositan 
# más arena y limo del glaciar en retroceso.

plot(astsa::varve)

logvarve = log(varve)

plot(logvarve)
acf2(logvarve,50)

#Contraste de Dickey-Fuller
adf.test(varve) #note que es engañoso.
adf.test(logvarve)
adf.test(diff(logvarve))

plot(diff(logvarve))
acf2(diff(logvarve),50)


astsa::sarima(logvarve, 0, 1, 1, no.constant=TRUE)   # ARIMA(0,1,1)

astsa::sarima(logvarve, 1, 1, 1, no.constant=TRUE)   # ARIMA(1,1,1)


