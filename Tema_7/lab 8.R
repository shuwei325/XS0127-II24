
# Laboratorio 8 -----------------------------------------------------------

library(ggplot2)
library(forecast)
library(fpp2)
library(astsa)
library(car)
library(TSA)
library(tseries)


# 1. Pronóstico del cambio de gasto basado en el ingreso personal  --------

?uschange

autoplot(uschange[,1:2], facets=TRUE) +
  xlab("Year") + ylab("") +
  ggtitle("Quarterly changes in US consumption
    and personal income")

y<-uschange[,"Consumption"]
x<-uschange[,"Income"]

#Ilustración de un modelo con errores ARIMA y un modelo de regresión en diferencias
mod0 <- Arima(y, xreg=x, order=c(1,1,0))
summary(mod0)

#equivalentes
mod1 <- Arima(diff(y), xreg=diff(x), order=c(1,0,0)) #asume que tiene intercepto
summary(mod1)
mod2 <- Arima(diff(y), xreg=diff(x), include.mean=FALSE , order=c(1,0,0))
summary(mod2)

#auto.arima
mod <- auto.arima(y, xreg=x)
summary(mod)

mod$residuals

checkresiduals(mod,lag=20)
res<-mod$res
ts.plot(res)
acf2(res)
shapiro.test(res)
qqPlot(res)

cbind("Regression Errors" = residuals(mod, type="regression"),
      "ARIMA errors" = residuals(mod, type="innovation")) %>%
  autoplot(facets=TRUE)

acf2(residuals(mod, type="regression"))
acf2(residuals(mod, type="innovation"))

#Estimación del modelo de regresión asumiendo errores independientes y normales
mod.reg<-lm(y~x)
summary(mod.reg)
res.lm<-mod.reg$residuals
ts.plot(res.lm)
acf2(res.lm)

mod.res <- auto.arima(res.lm)
summary(mod.res)

# 2. Tendencia determinística y estocástica -------------------------------

set.seed(123456)
e <- rnorm(500)
## caminata aleatoria
rw.nd <- cumsum(e)
## tendencia
tend <- 1:500
## caminata aleatoria con desvío
rw.wd <- 0.5*tend + cumsum(e)
## tendencia determinística con ruido
dt <- e + 0.5*tend
## plotting
par(mar=rep(5,4))
plot.ts(dt, lty=1, col=1, ylab='', xlab='')
lines(rw.wd, lty=2, col= 2)
par(new=T)
plot.ts(rw.nd, lty=3, col=3, axes=FALSE)
axis(4, pretty(range(rw.nd)))
legend(10, 18.7, legend=c('tend. determ. + ruido ',
                          'tend. determ. + tend. estocast.', 'tend. estocast.'),
       lty=c(1, 2, 3),col=c(1,2,3)) 


# 3. Regresión con variables independientes rezagadas ---------------------

#Cotización mensual y gastos en publicidad de una compañía estadounidense (enero, 2002- abril, 2005)

?insurance

autoplot(insurance, facets=TRUE) +
  xlab("año") + ylab("") +
  ggtitle("Cotización mensual y gastos en anuncios")

y<- insurance[,"Quotes"]
x<- insurance[,"TV.advert"]

# Predictores rezagadas (0,1,2 y 3 rezagos)

anuncios <- cbind(
  x0 = x,
  x1 = stats::lag(x,-1),
  x2 = stats::lag(x,-2),
  x3 = stats::lag(x,-3)) %>%
  head(NROW(insurance))                    #eliminar los NA al final


# Restringir datos para comparar los modelos del mismo periodo
(mod1 <- auto.arima(insurance[4:40,1], xreg=anuncios[4:40,1],
                    stationary=TRUE))
(mod2 <- auto.arima(insurance[4:40,1], xreg=anuncios[4:40,1:2],
                    stationary=TRUE))
(mod3 <- auto.arima(insurance[4:40,1], xreg=anuncios[4:40,1:3],
                    stationary=TRUE))
(mod4 <- auto.arima(insurance[4:40,1], xreg=anuncios[4:40,1:4],
                    stationary=TRUE))

c(mod1[["aicc"]],mod2[["aicc"]],mod3[["aicc"]],mod4[["aicc"]])


(mod.final <- auto.arima(insurance[,1], xreg=anuncios[,1:2],
                         stationary=TRUE))

#pronóstico

pronostico <- forecast(mod.final, h=20,
                       xreg=cbind(x0 = rep(8,20),
                                  x1 = c(anuncios[40,1], rep(8,19))))
autoplot(pronostico) + ylab("Cotización") +
  ggtitle("Proyección")



# Intervención ------------------------------------------------------------

# 4. Serie de colgate y crest --------

table<-read.table("crestcolgate.dat")
colnames(table)<-c( " CRESTMS", "COLGTEMS" ,"CRESTPR", "COLGTEPR")
names(table)

crest.colgate<-ts(table[,1:2])
autoplot(crest.colgate)

#Modelo ARIMA antes de la intervención
#las primeras 134 semanas:
dim(crest.colgate)
crest<-ts(table[,1])
crest1<-crest.colgate[1:134]
crest2<-crest.colgate[135:276]

ts.plot(crest)
points(135:276,crest2,col=2,type="l")
abline(v=134)
legend("topleft",c("antes","después"),lwd=2,col=1:2)


acf2(crest1,main="antes")

dif.crest1<-diff(crest1)
ts.plot(dif.crest1)
acf2(dif.crest1)

moda <- Arima(crest1, order=c(0,1,1))
summary(moda)

#creación de las variables indicadoras
I1 <- c(rep(0,134),rep(1,142)) #intervención durante la semana 135
I2 <- c(rep(0,135),rep(1,141)) #intervención durante la semana 136
X<-cbind(I1,I2)

modb <- Arima(crest, xreg=X, order=c(0,1,1))
summary(modb)


checkresiduals(modb,lag=20)
res<-modb$res
ts.plot(res)
acf2(res)

shapiro.test(res)
qqPlot(res)

