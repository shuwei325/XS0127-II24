
# Laboratorio 4 -----------------------------------------------------------

library(ggplot2)
library(forecast)
library(fpp2)
library(astsa)
library(car)

# 0. Modelo de Regresión --------------------------------------------------

?uschange
data(uschange)

autoplot(uschange[,1:5], facets = TRUE, colour=TRUE) +
  ylab("") + xlab("Year") +
  guides(colour="none")

uschange %>%
  as.data.frame() %>%
  GGally::ggpairs()

mod0 <- lm(
  Consumption ~ Income + Production + Unemployment + Savings,
  data=uschange)
summary(mod0)

mod0 <- tslm(
  Consumption ~ Income + Production + Unemployment + Savings,
  data=uschange)
summary(mod0)


autoplot(uschange[,'Consumption'], series="Data") +
  forecast::autolayer(fitted(mod0), series="Fitted") +
  xlab("Year") + ylab("") +
  ggtitle("Cambio porcentual de gastos de consumo en EE.UUU ") +
  guides(colour=guide_legend(title=" "))

checkresiduals(mod0)


# 1. Modelos de tendencia -------------------------------------------------

# 1.1. Ejemplo con graduados de ITCR de 1975 a 2002 -----------------------

itcrgrad<-read.csv("ITCR.csv",sep=",")
y<-ts(itcrgrad$graduados,start=1975)
autoplot(y) 

# crear variables independientes
tiempo<-seq(1,length(y))
tiempo2<-tiempo^2

datos.itcrgrad<-data.frame(y,tiempo,tiempo2)

mod1<-lm(y~tiempo+tiempo2,datos.itcrgrad)
summary(mod1)

mod.ts1<-tslm(y~trend+trend^2) #note la salida no incluye al tendencia al cuadrado.
summary(mod.ts1)
mod.ts2<-tslm(y~trend+I(trend^2))
summary(mod.ts2)

par(mfrow = c(1,1))
e<-mod1$residuals

#Normalidad
hist(e)
shapiro.test(e)

#Homoscedasticidad
ts.plot(e)
lmtest::bptest(mod1)

#Autocorrelación
lag1.plot(e)
acf(mod1$residuals)

durbinWatsonTest(mod1) #solamente considera un rezago

checkresiduals(mod1)

# Regresión con series estacionales ---------------------------------------

turistas<-read.csv("turistas.csv",sep=";")
y<-ts(turistas$turistas,start=c(1991,1),frequency=12)

#transformacion logarítmica
w<-log(y)
autoplot(y) 
autoplot(w)

tiempo<-seq(1,length(y))
tiempo2<-tiempo^2
mes<-rep(seq(1,12),10)
mes<-as.factor(mes)

datos1<-data.frame(w,tiempo,tiempo2,mes)
mod1<-lm(w~tiempo+tiempo2+mes,datos1)
summary(mod1)
levels(datos1$mes)

datos2 <- within(datos1, mes <- relevel(mes, ref = 12))
levels(datos2$mes)
mod2 <- lm(w~tiempo+tiempo2+mes,datos2)
summary(mod2)


#pronóstico
mod1<-lm(w~tiempo+tiempo2+mes,datos1)
summary(mod1)
pronostico<-forecast(mod1) #error

mod3<-tslm(w~trend+I(trend^2)+season)
summary(mod3)
pronostico<-forecast(mod3,h=12)
pronostico

autoplot(w) +
  ylab("ln Y") +
  autolayer(mod3$fitted.values, series = "ajustado") +
  autolayer(pronostico, series = "pronostico")

checkresiduals(mod3)

####

y<-ts(turistas$turistas,start=c(1991,1),frequency=12)
y.train<-window(y,start=c(1991,1),end=c(1999,12))
y.test<-window(y,start=c(2000,1),end=c(2000,12))
mod4<-tslm(y.train~trend+I(trend^2)+season)
summary(mod4)
pronostico<-forecast(mod4,h=12)
pronostico
accuracy(pronostico)
accuracy(pronostico,y.test)
