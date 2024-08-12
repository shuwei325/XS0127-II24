
# Laboratorio 7 -----------------------------------------------------------

library(ggplot2)
library(forecast)
library(fpp2)
library(astsa)
library(car)
library(TSA)
library(tseries)
library(sarima)
library(BETS)


# 0. Ejemplos simulados ---------------------------------------------------

#SARMA puro
TT<- 150 
# TT<- 500 (prueben con varios tamaños)

set.seed(1000)
#ej 1: SAR(1)_12
x1 <- sim_sarima(n=TT,model=list(sar=0.8, nseasons=12, sigma2 = 1))
ts.plot(x1)
acf2(x1,max.lag=50)

#las raíces del polinomio autoregresivo.
raices<-polyroot(c(1,rep(0,11), -0.8))
modulo.raices<-Mod(polyroot(c(1,rep(0,11), -0.8)))
modulo.raices>1

#ej 2: SMA(1)_4
x2 <- sim_sarima(n=TT,model=list(sma=0.8, nseasons=4, sigma2 = 1))  
ts.plot(x2)
acf2(x2,max.lag=30)

#ej 3: SARMA(1,1)_7
x3 <- sim_sarima(n=TT,model=list(sar=0.6,sma=0.8, nseasons=7, sigma2 = 1))  
ts.plot(x3)
acf2(x3,max.lag=30)


#SARMA mixto
TT<- 500 
# ej 1: SARMA(0,0,1)(1,0,0)_12
x1 <- sim_sarima(n=TT,model=list(ma=-.5,sar=0.8, nseasons=12, sigma2 = 1))
ts.plot(x1)
acf2(x1,max.lag=50)

# ej 2: SARMA(1,0,0)(0,0,1)_12
x2 <- sim_sarima(n=TT,model=list(ar=-.5,sma=0.8, nseasons=12, sigma2 = 1))
ts.plot(x2)
acf2(x2,max.lag=50)


# ej 3: SARMA(1,0,1)(1,0,0)_12
x3 <- sim_sarima(n=TT,model=list(ar=.7,ma=0.4,sar=0.9, nseasons=12, sigma2 = 1))
ts.plot(x3)
acf2(x3,max.lag=50)

# ej 3a: SARMA(12,0,1)(0,0,0)_12
x3a <- arima.sim(n = TT, model = list(order = c(12,0,0), 
                                      ar = c(0.7,rep(0,10),0.9)),sd=1,rand.gen= rnorm) + m
coef<-c(1,-0.7,rep(0,10),-0.9)
modulo.raices<-Mod(polyroot(coef))
modulo.raices>1
all(modulo.raices>1)

#devolviendo al ej 3: resolviendo (1-phi B)(1-Phi B^12)
coef<-c(1,-0.7,rep(0,10),-0.9,(-0.7*-0.9))
modulo.raices<-Mod(polyroot(coef))
modulo.raices>1
all(modulo.raices>1)
  
# ej 4: SARIMA(1,0,1)(0,0,1)_12
x4 <- sim_sarima(n=TT,model=list(ar=.7,ma=0.4,sma=0.9, nseasons=12, sigma2 = 1))
ts.plot(x4)
acf2(x4,max.lag=50)

# ej 5: SARIMA(1,0,1)(1,0,1)_12
x4 <- sim_sarima(n=TT,model=list(ar=.7,ma=0.4,sar=-0.5,sma=0.9, nseasons=12, sigma2 = 1))
ts.plot(x4)
acf2(x4,max.lag=50)

# ej 6: SARIMA(1,1,0)(0,1,1)_12
x5 <- sim_sarima(n=TT, model = list(ar=0.7,iorder=1, siorder=1,sma=0.4,
                                    nseasons=12, sigma2 = 1))
ts.plot(x5)
acf2(x5,max.lag=50)

x5<-x5[-c(1:14)]
ts.plot(x5)
acf2(x5,max.lag=50)

diff.x5 <- diff(x5,lag=1)
ts.plot(diff.x5)
acf2(diff.x5,max.lag=150)

sdiff.diff.x5 <- diff(diff.x5,lag=12)
acf2(sdiff.diff.x5,max.lag=50)


# 1. Serie mensual del número de turistas que ingresaron a Costa Rica --------

turistas<-read.csv("turistas.csv",sep=";")
y<-ts(turistas$turistas,start=c(1991,1),frequency=12)

#transformacion logarítmica
w<-log(y)
autoplot(y) 
autoplot(w)

acf2(w)

diffw<-diff(w)
autoplot(diffw)
acf2(diffw,max.lag=80)

sdif.diffw<-diff(diffw,lag=12)
autoplot(sdif.diffw)
acf2(sdif.diffw)

# ahora los correlogramas muestran estacionariedad.
# Identificación de los modelos apropiados.

# SARIMA(0,1,1)(0,1,1)_12
mod0a = Arima(sdif.diffw, order=c(0,0,1),
              seasonal=list(order=c(0,0,1),period=12))
summary(mod0a)
checkresiduals(mod0a,lag=20)
acf2(mod0a$res)
autoplot(mod0a)

mod0a$coef
coef.polinomio<-c(1,mod0a$coef[1],rep(0,10), mod0a$coef[2],mod0a$coef[1]*mod0a$coef[2])
raices<-polyroot(coef.polinomio)
modulo.raices<-Mod(raices)
modulo.raices>1
all(modulo.raices>1)

# SARIMA(1,1,0)(0,1,1)_12
mod0b = Arima(sdif.diffw, order=c(1,0,0),
              seasonal=list(order=c(0,0,1),period=12))
summary(mod0b)
checkresiduals(mod0b,lag=20)
acf2(mod0b$res)
autoplot(mod0b)


#Modelo final
mod.final = Arima(w, order=c(0,1,1),
                  seasonal=list(order=c(0,1,1),period=12))
summary(mod.final)
autoplot(mod.final)


res<-mod.final$res
ts.plot(res)
acf2(res)

#Normalidad
shapiro.test(res)
jarque.bera.test(res)
qqPlot(res)


# 2. Nivel de dióxido de carbono en alerta, Canada ------------------------
# Serie mensual del nivel de CO2 en alerta, de 01-1994 a 12-2004.

data(co2)
ts.plot(co2)
acf2(co2,50)

plot(diff(co2))
acf2(diff(co2),50)

plot(diff(diff(co2),lag=12))
acf2(diff(diff(co2,lag=12)),50)

#¿Cuáles son los modelos posibles?

mod1 = Arima(co2, order=c(0,1,1),
             seasonal=list(order=c(0,1,1),period=12))
summary(mod1)

t_test(model = mod1)

autoplot(mod1)
checkresiduals(mod1,lag=30)
tsdiag(mod1,gof.lag=30)


res<-mod1$res
ts.plot(res)
acf2(res)

#Normalidad
shapiro.test(res)
jarque.bera.test(res)

qqPlot(res)

astsa::sarima(co2, p=0,d=1,q=1, P=0, D=1, Q=1, S=12, no.constant=TRUE) 


# 3. Pasajeros de avión --------------------------

data(AirPassengers)
y <- AirPassengers

autoplot(y)

w<-log(y)
autoplot(w)

acf2(w)

diffw<-diff(w)
autoplot(diffw)
acf2(diffw)
acf2(diffw,max.lag=70)
acf2(diffw,max.lag=80)
#diferencia estacional

#vamos a ajustar sin diferencia estacional para ver qué pasa...
mod0a <- Arima(log(y), c(1, 1, 1),seasonal = list(order = c(1, 0, 0), period = 12))
summary(mod0a)
autoplot(mod0a)
checkresiduals(mod0a,lag=20)

mod0a <- Arima(log(y), c(1, 1, 1),seasonal = list(order = c(1, 0, 0), period = 12),method="CSS")
mod0a$coef
coef.polinomio<-c(1,-mod0a$coef[1],rep(0,10), -mod0a$coef[3],mod0a$coef[1]*mod0a$coef[3])
raices<-polyroot(coef.polinomio)
modulo.raices<-Mod(raices)
modulo.raices>1
all(modulo.raices>1)

acf2(mod0a$res)

#Diferencia estacional
sdif.diffw<-diff(diffw,lag=12)
autoplot(sdif.diffw)
acf2(sdif.diffw)

# ahora los correlogramas muestran estacionariedad.
# Identificación de los modelos apropiados.

# 1. Empezamos con la parte estacional
# SARIMA(0,1,0)(0,1,1)_12
mod1a <- Arima(log(y), c(0, 1, 0),seasonal = list(order = c(0, 1, 1), period = 12))
summary(mod1a)
checkresiduals(mod1a,lag=20)
acf2(mod1a$res)
autoplot(mod1a)

# SARIMA(0,1,0)(1,1,0)_12
mod1b <- Arima(log(y), c(0, 1, 0),seasonal = list(order = c(1, 1, 0), period = 12))
summary(mod1b)
checkresiduals(mod1b,lag=20)
acf2(mod1b$res)
autoplot(mod1b)

# 2. Luego seguimos con la parte no estacional

# SARIMA(1,1,0)(0,1,1)_12
mod2a <- Arima(log(AirPassengers), c(1, 1, 0),seasonal = list(order = c(0, 1, 1), period = 12))
summary(mod2a)
checkresiduals(mod2a,lag=24)
acf2(mod2a$res)
autoplot(mod2a)

# SARIMA(0,1,1)(0,1,1)_12
mod2b <- Arima(log(AirPassengers), c(0, 1, 1),seasonal = list(order = c(0, 1, 1), period = 12))
summary(mod2b)
checkresiduals(mod2b,lag=20)
acf2(mod2b$res)
autoplot(mod2b)



#Modelo final

mod.final <- mod2b
summary(mod.final)
autoplot(mod.final)

res<-mod.final$res
ts.plot(res)
acf2(res)

#Normalidad
shapiro.test(res)
jarque.bera.test(res)
qqPlot(res)

#Pronóstico
forecast(mod.final)
autoplot(forecast(mod.final))


#

# SARIMA(0,1,1)(0,1,1)_12
mod2b <- Arima(AirPassengers, c(0, 1, 1),seasonal = list(order = c(0, 1, 1), period = 12),lambda=0)
summary(mod2b)
autoplot(forecast(mod2b,h=50))

