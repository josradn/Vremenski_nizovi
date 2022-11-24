library(tseries)
library(TSA)
library(stats)
library(aTSA)
library(fGarch)
library(ggplot2)

setwd("C:/Users/Josipa/Pictures/Faks/DIPLOMSKI/VN/praktični")

data=read.table("VremenNiz-25los.txt")
data=t(data)
data=as.vector(data)
n=length(data)

#1)
ts.plot(data)
#Iz grafa možemo vidjeti da postoji trend i sezonalnost
ts.plot(data[20:50],gpars=list(ylab="Data"))
points(1:31,data[20:50],col="red",lwd=2)

#Pretpostavljam da postoji nekakva sezonalnost, koja bi mogla bit perioda 5
pod=log(data[2:n]/data[1:(n-1)])
ts.plot(pod,gpars=list(ylab="Data"))

#2)
#podaci u obliku timeseries
data_ts=ts(data,freq=4)
#podaci secirani po trendu sezonalnosti i sumu
data_stl=stl(data_ts,"periodic")
plot(data_stl)
#maknimo sezonalnost
data_s=data_stl$time.series[,1]
data_ds=data_ts-data_s

plot(data_ds)
ts.plot(data)
#acf(data)
#Vidimo da kad smo makli sezonalnost, da se podaci nisu poboljsali, dakle vjerojatno uopce nema sezonalnosti


#procjenimo trend log povratima
pod=log(data[2:n]/data[1:(n-1)])
ts.plot(pod,gpars=list(ylab="Data"))#vidimo da je uklonjen trend i nema sezonalnosti
#acf(pod)
#podaci izgledaju u redu, pa cemo raditi s tim novim podacima

#procijenimo trend polinomijalno (max 3 stupnja)
t=1:260
trend=lm(data~(t))
summary(trend)
trend2=lm(data~(t)+I(t^2))
summary(trend2)
anova(trend,trend2)#znacajno pa nam je trend2 bolji
trend3=lm(data~(t)+I(t^2)+I(t^3))
anova(trend2,trend3)#znacajno pa nam je trend3 bolji

ts.plot(data,gpars=list(ylab="Data"))
lines(t, trend$fitted.values,col="red",lwd=2)
lines(t, trend2$fitted.values,col="blue",lwd=2)
lines(t, trend3$fitted.values,col="green",lwd=2)
#Vidimo da polinom 3.stupnja najljepse fita trend

ts.plot(trend3$residuals,gpars=list(ylab="Data"))
data_dt=data-trend3$fitted.values#ili maknemo trend, ali to je isto sto i trend3$residuals
data_dt=trend3$residuals # vidimo da su podaci uglavnom rasprseni oko x osi i da je maknut trend i nema sezonalnosti


par(mfrow=c(1,2))
ts.plot(pod,gpars=list(ylab="Data",main="log-povrati"))
ts.plot(trend3$residuals,gpars=list(ylab="Data",main="polinom 3.stupnja"))

#3)
#Mozemo koristit data_dt ili pod s log povratima, rpvo cu napravit za data_dt
acf(data_dt)
#Vidimo iz grafa da bi mogli naslutiti da se radi o MA(4) ili MA(8) ili MA(9), no kako necemo radit veci od 3 vjerojatno ce on bit najbolji
pacf(data_dt)
#Vidimo iz grafa da bi najbolji bio AR(3)

#Tesko je zakljucivati iz acf i pacf za modele, pa cemo za AR(p) model zakljuciti Yule-Walker metodom
ar_p=ar(data_dt,method="yule-walker")
summary(ar_p)
p=ar_p$order
p
ar_p_koef=ar_p$ar
#Vidimo iz Yule-Walker metode da je najbolje uzeti model AR(4)

#Sada cemo procijeniti parametar za MA(q) model pomocu AIC kriterija
ma_1=arima(data_dt,order=c(0,0,1),include.mean = F)
ma_1$aic #2713.776

ma_2=arima(data_dt,order=c(0,0,2),include.mean = F)
ma_2$aic #2712.533

ma_3=arima(data_dt,order=c(0,0,3),include.mean = F)
ma_3$aic #2702.016

ma_4=arima(data_dt,order=c(0,0,4),include.mean = F)
ma_4$aic #2700.663
#Vidimo pomocu AIC kriterija da je od MA modela najbolji s parametrom 3, tj MA(3), sto smo i pretpostavili iz grafa, cak bi MA(4) bio bolji al necemo to gledat


#Sada cemo izabrati najbolji ARMA(p,q) model pomocu AIC kriterija
arma_AIC <- function(vr_niz,n){
  min = arima(vr_niz,order=c(0,0,0),include.mean = T)$aic
  p1 = 0
  q1 = 0
  for (i in 0:n){
    for (j in 0:(n-i)){
      tren = arima(vr_niz,order=c(i,0,j),include.mean = T)$aic
      if (tren < min){
        min = tren
        p1 = i
        q1 = j
      }
    }
  }
  return(c(p1,q1,min))
}
arma_AIC(data_dt,7)
#ARMA(3,3),2697.657
arma_AIC(data_dt,6)
#ARMA(3,3),2697.657
arma_AIC(data_dt,5)
#ARMA(1,2),2699.635
arma_AIC(data_dt,4)
#ARMA(1,2),2699.635
#Mozemo cak reci po AIC kriteriju da bi bolji bio model ARMA(3,3) no vidit cu ocu taj zakljucak ukljucit


#Pomocu AIC kriterija cemo sada usporediti modele AR(4), ARMA(1,1) i MA(3)
ar_4=arima(data_dt,order=c(4,0,0),include.mean = F)
ar_4$aic #2698.197
ma_3$aic #2702.016
arma_11=arima(data_dt,order=c(1,0,1),include.mean = F)
arma_11$aic #2698.277
#Pomocu AIC kriterija zakljucujemo da je od ta 3 modela najbolji MA(3) model.
#Takoder ako ga usporedimo s nekim od gore ARMA modela, takoder je bolji po AIC kriteriju od svih.
#Dakle, po AIC kriteriju MA(3) je najbolji model. 


#Sada cu iste zakljucke izvest i za podatke s logpovratima
#mozemo koristiti pod ili data_dt, meni se cini da je bolje s log povratima pa cu koristit za njih (mozda napravim oboje)
acf(pod)
#Vidimo iz acf grafa da bi mogao biti model MA(1)
pacf(pod)
#Vidimo iz pacf grafa da bi mogao biti model AR(4) ili AR(7)

#Tesko je zakljucivati iz acf i pacf za modele, pa cemo za AR(p) model zakljuciti Yule-Walker metodom
ar_pl=ar(pod,method="yule-walker")
summary(ar_pl)
pl=ar_pl$order
pl
ar_pl_koef=ar_pl$ar
#Vidimo iz Yule-Walker metode da je najbolje uzeti model AR(7), sto potvrduje predvidanja iz pacf modela


#Sada cemo procijeniti parametar za MA(q) model pomocu AIC kriterija
ma_1l=arima(pod,order=c(0,0,1),include.mean = F)
ma_1l$aic #555.8641

ma_2l=arima(pod,order=c(0,0,2),include.mean = F)
ma_2l$aic #556.9661

ma_3l=arima(pod,order=c(0,0,3),include.mean = F)
ma_3l$aic #558.9158
#Vidimo pomocu AIC kriterija da je model MA(1) najbolji, sto potvrduje pretpostavke iz acf grafa


#Sada cemo izabrati najbolji ARMA(p,q) model pomocu AIC kriterija
arma_AIC(pod,10)
#ARMA(3,5),553.4992
arma_AIC(pod,9)
#ARMA(3,5),553.4992
arma_AIC(pod,8)
#ARMA(3,5),553.4992
arma_AIC(pod,7)
#ARMA(5,1),554.3738
arma_AIC(pod,6)
#ARMA(5,1),554.3738
#Prema AIC kriteriju najbolji bi bio ARMA(3,5) model, no jos cemo vidit hocemo li njega ukljucit

#Pomocu AIC kriterija cemo sada usporediti modele AR(7), ARMA(1,1) i MA(1)
ar_7l=arima(pod,order=c(7,0,0),include.mean = F)
ar_7l$aic #568.545
ma_1l$aic #555.8641
arma_11l=arima(pod,order=c(1,0,1),include.mean = F)
arma_11l$aic #556.954
#Pomocu AIC kriterija zakljucujemo da je od ta 3 modela najbolji MA(1) model.
#Ako ga usporedimo nekim od gore navedenih ARMA modela, ipak je ARMA(3,5) bolji, no mozda me to ne pita.
#Dakle, po AIC kriteriju MA(1) je najbolji model od ta 3 navedena.




#4)
#Pod 3) sam vec odgovorila zasto biram MA(1) model za logpovrate, i MA(3) model za polinomijalno
#Sada cu isprintati koeficijente navedenih modela:
#polinomijalno -> Xt=Zt+0.04511168*Zt-1+0.10900138*Zt-2+0.23063933*Zt-3, sigma2=1864.406
ma_3$coef
ma_3$sigma2
acf(ma_3$residuals)#trebali bi bit svi unutar pruge, sto se vidi da su 2/24 samo vani pa mozemo zakljuciti da je model dobar
#reziduali=fitted-stvarni -> fitted=stvarni+reziduali
ma_3fit=data_dt+ma_3$residuals
ts.plot(data_dt)
ts.plot(ma_3fit)
Box.test(ma_3$residuals, type="Ljung") # testira jesu li podaci IID(0,sigma^2)
#dobila sam pvrj=0.6526, dakle ne odbacujemo hipotezu da podaci jesu IID(0,sigma^2), dakle bijeli sum su WN(0,sigma^2)


#logpovrati -> Xt=Zt-0.8889*Zt-1, sigma2=0.493894
ma_1l
ma_1l$coef
ma_1l$sigma2
acf(ma_1l$residuals)#navodno svi trebaju biti unutar pruge, sto se vidi i da jesu osim 1/24 pa mozemo zakljuciti da je model dobar
ma_1lfit=pod+ma_1l$residuals
ts.plot(ma_1lfit)
ts.plot(pod)
Box.test(ma_1l$residuals, type="Ljung")
#dobila sam pvrj=0.3981, dakle ne odbacujemo hipotezu da podaci jesu IID(0,sigma^2), dakle bijeli sum su WN(0,sigma^2)


#3) JOS TREBA RASPRAVITI GARCH(1,1)
#polinomijalno:
acf(abs(pod))
acf(pod^2)

arch.test(ma_3,output=TRUE)
#sve su p vrijednosti niske, pa odbacujemo hipotezu H0 o homoskedasticnosti
#sve su p vrijednosti visoke, pa opet potvrdujemo da su reziduali bijeli sum

#logpovrati:
arch.test(ma_1l,output=TRUE)
#sve su p vrijednosti 0, pa odbacujemo hipotezu H0 o homoskedasticnosti
#sve su p vrijednosti visoke, pa opet potvrdujemo da su reziduali bijeli sum 
McLeod.Li.test(ma_1l,pod)
#Dakle, GARCH(1,1) model cemo sad fitat jer ocito su podaci heteroskedasticni, vidjeti cemo kako on dobro fita podatke:
garchl=garchFit(formula=~garch(1,1),data=pod)
garchlfit=garchl@residuals+pod
ts.plot(garchlfit)
ts.plot(pod)

garchp=garchFit(formula=~garch(1,1),data=data_dt)
garchpfit=garchp@residuals+data_dt
ts.plot(garchpfit)
ts.plot(data_dt)


#5)

#Za polinomijalno koristit cemo model ma_3, a za log povrate koristit cemo ma_1l

#Ovo dobivamo uz pretpostavku uz to da je niz gaussovski (stac s očekivanjem 0)??
#polinomijalno:
new_vrjp=predict(ma_3, n.ahead=1)$pred[1]
new_vrjp

Up = new_vrjp + 1.96*predict(ma_3, n.ahead=1)$se[1]
Lp = new_vrjp - 1.96*predict(ma_3, n.ahead=1)$se[1]
Up
Lp

#logpovrati
new_vrjl=predict(ma_1l, n.ahead=1)$pred[1]
new_vrjl
prava_new=exp(new_vrjl)*data[260]
prava_new

Ul = new_vrjl + 1.96*predict(ma_1l, n.ahead=1)$se[1]
Ll = new_vrjl - 1.96*predict(ma_1l, n.ahead=1)$se[1]
Ul
Ll
prava_ul=exp(Ul)*data[260]
prava_ll=exp(Ll)*data[260]
prava_ul
prava_ll

forecast(ma_3,lead=1,alpha=0.05)# isto je JEJ
forecast(ma_1l,lead=1,alpha=0.05)# isto je JEJ


par(mfrow=c(1,2))
df=data.frame(c(data[251:260],prava_new))
df
x_predikcija=251:261
y_predikcija=df$c.data.251.260...prava_new.
y_donji=c(y_predikcija[1:10],prava_ll)
y_gornji=c(y_predikcija[1:10],prava_ul)
ggplot(df,aes(x=x_predikcija)) + 
  geom_line(aes(y=y_predikcija,col="Stvarni podaci + predikcija"),colour='blue',size=1.2,show.legend = TRUE)  + 
  geom_line(aes(y=y_donji,col="Donja granica 95% pouzdanog intervala"),colour='red',linetype = "dashed")+ 
  geom_line(aes(y=y_gornji,col="Gornja granica 95% pouzdanog intervala"),colour='red',linetype = "dashed")+
  ggtitle("Predikcija i 95% pouzdani interval")



plot(251:260,data[251:260],type="l",xlim=c(251,261),ylim=c(0,600),xlab="Time",ylab="Data")
#plot(701:710,data.bez.trend[701:710],type="l",xlim=c(701,711))
lines(261,prava_new, col="red", type="b",lwd=2)
lines(261,prava_ul, col="blue",type="o",lwd=2)
lines(261,prava_ll, col="blue",type="o",lwd=2)
