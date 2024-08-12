library(gamair)
library(mgcv)
data(cairo)

ctamm <- gamm(temp~s(day.of.year,bs="cc",k=20)+s(time,bs="cr"),
              data=cairo,correlation=corAR1(form=~1|year))

pred<-fitted(ctamm$gam)
plot(cairo$temp~cairo$time,type="l", ylab="temp",xlab="time")
points(cairo$time,pred,type="l",col=2,lwd=2)


plot(cairo$temp,xlab='Day',ylab='Brightness')

library(astsa)
x<-cairo$temp
x.spec <-mvspec(x,log="no")
plot(x.spec)
frecuencia<-x.spec$freq[x.spec$spec==max(x.spec$spec)]
abline(v=frecuencia,col=2)
1/frecuencia

cos1<-cos(2*pi*frecuencia*cairo$time)
sin1<-sin(2*pi*frecuencia*cairo$time)
data=data.frame(x=x,time=cairo$time,cos=cos1,sen=sin1)
mod<-lm(x~cos1+sin1,data=data)
plot(cairo$temp~cairo$time,type="l", ylab="temp",xlab="time")
points(cairo$time,fitted(mod),type="l",col=2,lwd=2)

