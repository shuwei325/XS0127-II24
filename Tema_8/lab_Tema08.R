
# Tema VIII: Modelos lineales multivariados de series temporales ---------------

library(ggplot2)
library(forecast)
library(fpp2)
library(astsa)
library(car)
library(TSA)
library(tseries)
library(sarima)
library(vars)
library(dse) 

# 0. Ejemplos (descriptivo) ---------------------------------------------------


par(mfrow = c(2,1))
tsplot(soi, ylab="", main="Southern Oscillation Index")
tsplot(rec, ylab="", main="Recruitment") 


par(mfrow=c(2,1))
acf1(soi, 48, main="Southern Oscillation Index")
acf1(rec, 48, main="Recruitment")


lag2.plot (soi, rec, 10) #importante ver la dirección que está calculando
lag2.plot (soi, rec, 5)

par(mfrow=c(2,1))
ccf2(rec,soi, 48, main="función de correlación cruzada de Rec. contra SOI")
ccf2(soi,rec, 48, main="función de correlación cruzada de SOI contra Rec.")

(r1=ccf(rec,soi, 10, plot=FALSE))
(r2=ccf(soi,rec, 10, plot=FALSE))



# 1. Ejemplo de VAR(2) simulado ----------------------------------------------
#Code 2.1 tomado del Bernhard 2008

## Especificar los parámetros del polinomio autoregresivo: A(L) 
Apoly   <- array(c(1.0, -0.5, 0.3, 
                   0, 0.2, 0.1, 
                   0, -0.2, 0.7, 
                   1, 0.5, -0.3) ,
                 c(3, 2, 2))
#I+Phi_1 B + Phi_2 B
Apoly[1,,]
Apoly[2,,]
Apoly[3,,]
## Especificar la estructura de las inovaciones.
B <- diag(2)
## Especificar el término constante
TRD <- c(5, 10)
## Generación de VAR(2) 
var2  <- ARMA(A = Apoly, B = B, TREND = TRD)
## Simulación de 500 observaciones
varsim <- simulate(var2, sampleT = 500,
                   noise = list(w = matrix(rnorm(1000),
                                           nrow = 500, ncol = 2)), rng = list(seed = c(123456))) 
## Obtención de la serie generada
vardat <- matrix(varsim$output, nrow = 500, ncol = 2)
colnames(vardat) <- c("y1", "y2")

plot.ts(vardat, main = "", xlab = "")

acf(vardat)

## Determinar el lag de acuerdo a los criterios
infocrit <- VARselect(vardat, lag.max = 3,
                      type = "const")
infocrit

## Estimar el modelo
varsimest <- VAR(vardat, p = 2, type = "const",
                 season = NULL, exogen = NULL)

summary(varsimest)

## Selección de acuerdo a BIC
varsimest <- VAR(vardat, type = "const",
                 lag.max = 3, ic = "SC")
## Verificar que los eigenvalues tengan módulo menor a 1 (estacionariedad)
roots <- vars::roots(varsimest)
roots

# 2. Ejemplo: contaminación, temperatura y mortalidad ------------------------
#Ejemplo 5.10 de Shumway & Stoffer
?cmort
?tempr
?part

par(mfrow=c(3,1))
tsplot(cmort, main="Cardiovascular Mortality", ylab="")
tsplot(tempr, main="Temperature",  ylab="")
tsplot(part, main="Particulates", ylab="")

pairs(cbind(Mortality=cmort, Temperature=tempr, Particulates=part))

acf2(cmort, 200, main="")
acf2(tempr, 100, main="")
acf2(part, 100, main="")

par(mfrow=c(3,1))
ccf2(cmort,tempr, 50, main="cmort vs tempr")
ccf2(cmort,part, 50, main="cmort vs part")
ccf2(tempr,part, 50, main="tempr vs part")

data = data.frame(cmort, tempr, part)
plot.ts(data , main = "", xlab = "")

stats::acf(data,lag.max = 100)

#VAR(1)
modvar1= VAR(data, p=1, type="both")
summary(modvar1)

acf(residuals(modvar1)[,1])
acf(residuals(modvar1)[,2])
acf(residuals(modvar1)[,3])

#VAR(2)

modvar2 = VAR(data, p=2, type="both")
summary(modvar2)

acf(residuals(modvar2)[,1],lag.max=100)
acf(residuals(modvar2)[,2])
acf(residuals(modvar2)[,3])

stats::acf(residuals(modvar2))

VARselect(data, lag.max=10, type="both")

#AIC (p=9)
#HQ=Hannan-Quinn (p=5)
#SC=BIC o criterio de información de Schwarz (p=2)  
#FPE=Final Predictor Error (p=9)

mod.f1 <- VAR(data, p=2, type="both")
mod.f2 <- VAR(data, p=5, type="both")
mod.f3 <- VAR(data, p=9, type="both")

summary(mod.f1)
summary(mod.f2)
summary(mod.f3)
stats::acf(resid(mod.f1), 52)
stats::acf(resid(mod.f2), 52)
stats::acf(resid(mod.f3), 52)

#Versión multivariada del contraste de Ljung-Box
serial.test(mod.f1, lags.pt=12, type="PT.adjusted")
serial.test(mod.f2, lags.pt=12, type="PT.adjusted")
serial.test(mod.f3, lags.pt=12, type="PT.adjusted")

#Predicción
(mod.pronostico = predict(mod.f1, n.ahead = 24, ci = 0.95))  # 4 semanas

fanchart(mod.pronostico)

# 3. Ejemplo: crecimiento de Producto Interno Bruto de UK, Canada y US ------------------------
#Ejemplo 4.7 Tsay


data=read.table("q-gdp-ukcaus.txt",header=T)
names(data)
pib=log(data[,3:5])
plot.ts(pib , main = "", xlab = "")

library(MTS)
z=diffM(pib)
z=z*100
stats::acf(z,lag.max = 100)


VARselect(z, lag.max=10, type="both")

modvar <- vars::VAR(z, p=2, type="both")

summary(modvar)
stats::acf(resid(modvar), 52)

serial.test(modvar, lags.pt=12, type="PT.adjusted")

shapiro.test(residuals(modvar)[,1])
shapiro.test(residuals(modvar)[,2])
shapiro.test(residuals(modvar)[,3])
mvnormtest::mshapiro.test(t(residuals(modvar)))

(mod.pronostico = predict(modvar, n.ahead = 5, ci = 0.95))  

fanchart(mod.pronostico)


