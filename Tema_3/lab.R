
# Laboratorio 3 -----------------------------------------------------------

library(ggplot2)
library(forecast)
library(fpp2)
library(plotly)


# Métodos de suavizamiento exponencial ------------------------------------


# 1. Suavizamiento exponencial simple ---------------------------------------

defunciones<-read.csv("defunciones.csv",sep=",")
y<-ts(defunciones$defunciones,start=c(2001,1),frequency=12)
autoplot(y) 
ses1 <- ses(y,alpha=0.3,initial="simple")
names(ses1)
ses1$model
ses1$fitted #valores ajustados
ses1$residuals #residuos


#estimación de alpha

alpha<-seq(0,1,0.01)
rmse<-as.numeric()
for(i in 1:length(alpha)){
  ses1 <- ses(y,alpha=alpha[i],initial="simple")
  rmse[i]<-accuracy(ses1)[2]
}
rmse
plot(alpha,rmse,type = "l")
alpha[rmse==min(rmse)]
abline(v=alpha[rmse==min(rmse)])



ses2 <- ses(y,h=12,initial="simple")
names(ses2)
ses2$model
ses2

autoplot(ses2)

autoplot(ses2) +
  autolayer(fitted(ses2), series="ajustado") +
  ylab("defunciones") + xlab("mes")


# 2. Método de Holt -------------------------------------------------------


itcrgrad<-read.csv("ITCR.csv",sep=",")
y<-ts(itcrgrad$graduados,start=1975)
autoplot(y) 
holt1 <- holt(y,alpha=0.3 ,beta=0.64 , h=5, initial="simple")
holt1

names(holt1)
holt1$model
holt1$mean # pronóstico con 5 pasos para adelante.
holt1$x # serie original
holt1$fitted # valores ajustados


holt2 <- holt(y , h=5, initial="simple")
holt2
names(holt2)
holt2$model

cbind(y,holt2$fitted)

autoplot(holt2) +
  autolayer(fitted(holt2), series="ajustado") +
  ylab("defunciones") + xlab("mes")


#estimación de alpha y beta

alpha<-seq(0,1,0.01)
beta<-seq(0,1,0.01)
rmse<-matrix(NA,nrow=length(alpha),ncol=length(beta))
for(i in 1:length(alpha)){
  for(j in 1:length(beta)){
    ses1 <- holt(y,alpha=alpha[i],beta=beta[j],initial="simple")
    rmse[i,j]<-accuracy(ses1)[2]
  }
}
rmse
rownames(rmse)<-alpha
colnames(rmse)<-beta


p <- plot_ly(z = rmse,x=alpha,y=beta, type = "surface")
p %>% layout(scene = list(xaxis = list(title = 'alpha'),
                          yaxis = list(title = 'beta'),
                          zaxis = list(title = 'rmse')))
#Debe abrir el gráfico en una ventana por a parte.

par.min<-which(rmse==min(rmse),arr.ind=TRUE)
alpha.min<-alpha[par.min[1]]
beta.min<-beta[par.min[2]]
holt1<-holt(y,alpha=alpha.min,beta=beta.min , initial="simple")
holt1$model
accuracy(holt1)

holt2<-holt(y , initial="simple")
holt2$model
accuracy(holt2)

holt3<-holt(y , initial="optimal")
holt3$model
accuracy(holt3)



# Comparación Holt - Holt amortiguado
fc <- holt(y, h=15)
fc2 <- holt(y, damped=TRUE, phi = 0.9, h=15)
fc2 <- holt(y, damped=TRUE, h=15)
ggplot2::autoplot(y) +
  autolayer(fc, series="Holt", PI=FALSE) +
  autolayer(fc2, series="Holt amortiguado", PI=FALSE) +
  ggtitle('Pronóstico usando método Holt-Holt amortiguado') + xlab("Año") +
  ylab("Graduados del ITCR (1975-2002)") +
  guides(colour=guide_legend(title="Pronóstico"))


# 3. Método multiplicativo de Holt-Winters ---------------------------------

turistas<-read.csv("turistas.csv",sep=";")
y<-ts(turistas$turistas,start=c(1991,1),frequency=12)
ts.plot(y)

ht1 <- hw(y,seasonal="multiplicative")
names(ht1)
ht1$model

ht2 <- hw(y,seasonal="additive")
ht3 <- hw(y, damped=TRUE, seasonal="multiplicative")

ht3$model

autoplot(y) +
  autolayer(ht1, series="HW multiplicativo", PI=FALSE, size = 1.2) +
  autolayer(ht2, series="HW aditivo", PI=FALSE, size = 1.2) +
  autolayer(ht3, series="HW M. amortiguado", PI=FALSE, size = 1.2) +
  xlab("fecha") +
  ylab("turistas") +
  ggtitle("Turistas que ingresaron a CR") +
  guides(colour=guide_legend(title="pronóstico"))

accuracy(ht1)
accuracy(ht2)
accuracy(ht3)



# 4. ETS como modelo Espacio de Estados -----------------------------------


fit4 <- ets(y)
summary(fit4)
autoplot(fit4)

plot(forecast(fit4))

fit4 %>% forecast(h=12) %>%
  autoplot() +
  ylab("Turistas")

