source("data.R")
source("functions.R")
library(KernSmooth)




#_____________________________________________________________________________________#

#figure 1
par(mfrow = c(2, 2)) 
#f1 haut a gauche
n=length(data1)
x=seq(0,n,1)[1:n]
plot(x/n,data1,pch=19,cex=.5,
     main= 'f1 : n = 100 et sigma = 0.1',
     xlab='x',ylab='f1(x)')
curve(f1,lwd=2,add=TRUE)
#f2 haut a droite
n=length(data2)
x=seq(0,n,1)[1:n]
plot(x/n,data2,pch=19,cex=.5,
     main= 'f2 : n = 100 et sigma = 0.1',
     xlab='x',ylab='f2(x)')
curve(f2,lwd=2,add=TRUE)
#f3 bas a gauche
n=length(data3)
x=seq(0,n,1)[1:n]
plot(x/n,data3,pch=19,cex=.5,
     main= 'f3 : n = 250 et sigma = 1',
     xlab='x',ylab='f3(x)')
curve(f3,lwd=2,add=TRUE)
#f4 bas a droite
n=length(data4)
x=seq(0,n,1)[1:n]
plot(x/n,data4,pch=19,cex=.5,
     main= 'f4 : n = 200 et sigma = 3',
     xlab='x',ylab='f4(x)')
curve(f4,lwd=2,add=TRUE)


#_____________________________________________________________________________________#

#figure 1.1
par(mfrow = c(2, 2)) 
#f1 haut a gauche
plot.error.function(f1,data1,0.1,"f1")
#f2 haut a droite
plot.error.function(f2,data2,0.1,"f2")
#f3 bas a gauche
plot.error.function(f3,data3,1,"f3")
#f4 bas a droite
plot.error.function(f4,data4,3,"f4")

#_____________________________________________________________________________________#

#figure 1.2
f.func.N(data1,f1,.1)

#_____________________________________________________________________________________#

#figure 1.3
f.func.N(data4,f4,3) 
#_____________________________________________________________________________________#


#figure 1.4
par(mfrow=c(2,2))
plot.sigma(data1,"f1")
plot.sigma(data2,"f2")
plot.sigma(data3,"f3")
plot.sigma(data4,"f4")


#_____________________________________________________________________________________#

#figure 2.1

h=0.1
n=100
t=seq(0,1,1/n)[1:n+1]
par(mfrow = c(1, 1)) 
alea=c(0.2,0.25,0.6,0.61,.59,1,.95,.9,.85,0)
exemple=matrix(0,n,10)

for(i in (1:10)){
  exemple[,i]=W.NW(t,n,alea[i]*100,K.gauss,h)
}
plot(t,apply(exemple,1,sum),col='brown',type='l',lwd=2,ylim=c(0,0.2),
     main="Fonction obtenue par somme de noyaux gaussiens", ylab="f(t)")
for(i in (1:10)){
  points(t,exemple[,i],type='l',lty=i,col='blue')
}



#_____________________________________________________________________________________#

#figure 2.2
par(mfrow = c(2, 2)) 
h=0.1
n=length(data4)
x=seq(0,1,1/n)[1:n+1]
plot(x,f4(x),type='l',lwd=2,lty=2,main='f1')
points(x,f.NW(data4,K.gauss,h),col=2,type='l',lwd=2)
legend("topleft",legend='f NW',col=2,lwd=2)

plot(x,f4(x),type='l',lwd=2,lty=2,main='f2')
points(x,f.NW(data4,K.rect,h),col=3,type='l',lwd=2)
legend("topleft",legend='f NW',col=3,lwd=2)

plot(x,f4(x),type='l',lwd=2,lty=2,main='f3')
points(x,f.NW(data4,K.trian,h),col=4,type='l',lwd=2)
legend("topleft",legend='f NW',col=4,lwd=2)

plot(x,f4(x),type='l',lwd=2,lty=2,main='f4')
points(x,f.NW(data4,K.parab,h),col=6,type='l',lwd=2)
legend("topleft",legend='f NW',col=6,lwd=2)
#_____________________________________________________________________________________#

#figure 2.3

par(mfrow=c(1,1))
n=1000
i=500
x=seq(0,1,1/n)
N=c(5,10,15,20)

plot(x,W_proj(n,i,x,N[4]),xlim=c(0,1),type='l',lwd=2,ylab="W_ni",main="Poids de projection")
for(k in (2:4)){
  points(x,W_proj(n,i,x,N[k-1]),type='l',col=c(2,4,6)[k-1],lwd=2,lty=k)
}
legend("topright",legend=c("N = 20","N = 15","N = 10","N = 5"),
       lty=c(1:4),lwd=2,col=c(1,2,4,6))


#_____________________________________________________________________________________#



#figure 2.4
par(mfrow = c(2, 2)) 
#graphe étude erreur data1
H1=seq(0,0.1,0.005)[-1]
Comp.fenetre(data1,H1,f1)
#graphe étude erreur data2
H2=seq(0,.1,0.01)[-1]
Comp.fenetre(data2,H2,f2)
#graphe étude erreur data3
H3=seq(0,.025,0.002)[-1]
Comp.fenetre(data3,H3,f3)
#graphe étude erreur data4
H4=seq(0,.05,0.005)[-1]
Comp.fenetre(data4,H4,f4)


#_____________________________________________________________________________________#

#figure 2.5
par(mfrow=c(2,2))
#f1 haut a gauche
Comp.estim(f1,data1,'f1')
#f2 haut a droite
Comp.estim(f2,data2,'f2')
#f3 bas a gauche
Comp.estim(f3,data3,'f3')
#f4 bas a droite
Comp.estim(f4,data4,'f4')



#_____________________________________________________________________________________#


#figure 3.1
par(mfrow = c(1, 2))
x=seq(0,1,1/length(data))[1:length(data)]
plot(data,pch=19,cex=.5,main='Moyenne mensuelle des températures',xlab="Mois",ylab="temperature")
plot.sigma(data,"Estimation de l'ecart-type des données ")

#_____________________________________________________________________________________#

#tableau du N mallows en fonction de sigma
sigma=c(1,2,3,4,8,9,10)
Nmallows=sapply( sigma , function(sigma){N.mallows(data, sigma)})
print(t(matrix(c(sigma,Nmallows),7,2)))


#_____________________________________________________________________________________#

#figure 3.2
par(mfrow = c(1, 2))
x=seq(0,1,1/length(data))[1:length(data)]
plot(data,pch=19,cex=.5,main='Estimateur pour s = 4 (N=42)',
     xlab="moment multipolaire",ylab="Spectre de puissance TT")
points(f.estim(data,42)(x),type='l',col='red',lwd=2)
plot(data,pch=19,cex=.5,main='Estimateur pour s = 2 (N=82)',
     xlab="moment multipolaire",ylab="Spectre de puissance TT")
points(f.estim(data,82)(x),type='l',col='red',lwd=2)
#erreur
par(mfrow = c(2, 2))
x=seq(0,1,1/length(data))[1:length(data)+1]
g1 = f.estim(data,42)(x)
g2 = f.estim(data,82)(x)
e1=data-g1
e2=data-g2  
plot(e1,type='l',main='vraies données : s = 4')
plot(e2,type='l',main='vraies données : s = 2')
plot(x=g1,y=e1,pch=19,cex=.5,main='vraies données : s = 4')
plot(x=g2,y=e2,pch=19,cex=.5,main='vraies données : s = 2')


#_____________________________________________________________________________________#

#figure 3.3
#centrer les donnees
databis=data-mean(data)
par(mfrow=c(1,1))
n=length(databis)
x=seq(0,1,1/n)[-(n+1)]
plot(x,databis,pch=19,cex=.5,col="darkgrey",
     main="moyennes mensuelles des températures",
     xlab="mois",ylab="temperature")
h=dpill(x,databis)
points(x,f.NW(databis,K.gauss,h),lwd=2,type='l',col=4)
legend("topleft",col=4,legend="estimateur NW pour h plug-in",lwd=2)
#erreur de NW
g = f.NW(databis,K.gauss,h)
e=databis-g
par(mfrow = c(1, 2))  
plot(e,type='l')
plot(x=g,y=e,pch=19,cex=.5)


#_____________________________________________________________________________________#

#figure 3.4
par(mfrow=c(1,2))
plot(datawave,pch=19,cex=.5,main="WMAP première année",
     xlab="moment multipolaire",ylab="Spectre de puissance TT")
plot.sigma(datawave,"Estimation de l'ecart-type des données ")


#_____________________________________________________________________________________#

#tableau associe
sigma=seq(900,1200,50)
Nmallows=sapply( sigma , function(sigma){N.mallows(datawave, sigma)})
print(t(matrix(c(sigma,Nmallows),7,2)))


#_____________________________________________________________________________________#

#figure 3.5
par(mfrow=c(1,2))
n=length(datawave)
x=seq(0,1,1/n)[-(n+1)]
plot(x,datawave,pch=19,cex=.5,col="darkgrey",
     main="Estimateur pour s = 1000 (N=14)",
     xlab="moment multipolaire",
     ylab="Spectre de puissance TT")
points(x,f.estim(datawave,14)(x),lwd=2,type='l',col=4)
plot(x,datawave,pch=19,cex=.5,col="darkgrey",
     main="Estimateur pour s = 1100 (N=7)",
     xlab="moment multipolaire",
     ylab="Spectre de puissance TT")
points(x,f.estim(datawave,7)(x),lwd=2,type='l',col=4)
#etude de l'erreur
par(mfrow = c(2, 2)) 
  g1 = f.estim(datawave,14)(x)
  g2 = f.estim(datawave,7)(x)
  e1=datawave-g1 
  e2=datawave-g2 
  plot(e1,type='l',main="vraies donnees : s = 1000")
  plot(e2,type='l',main="vraies donnees : s = 1100")
  plot(x=g1,y=e1,pch=19,cex=.5,main="vraies donnees : s = 1000")
  plot(x=g2,y=e2,pch=19,cex=.5,main="vraies donnees : s = 1100")
  

  
#_____________________________________________________________________________________#

#figure 3.6
par(mfrow=c(1,1))
n=length(datawave)
x=seq(0,1,1/n)[-(n+1)]
  plot(x,datawave,pch=19,cex=.5,col="darkgrey",
       main="Comparaison des deux estimateurs",
       xlab="moment multipolaire",ylab="Spectre de puissance TT")
  hwave=dpill(x,datawave)
  points(x,f.NW(datawave,K.gauss,hwave),lwd=2,type='l',col=2)
  points(x,f.estim(datawave,7)(x),lwd=2,type='l',col=4)
  legend("topleft",col=c(2,4),lwd=2,legend=c("estimateur NW pour h plug-in",
                                       "estimateur par projection pour N = 7"))

  
#_____________________________________________________________________________________#
  
  #annexes :
#estimateur projection pour differents N de f2
f.func.N(data2,f2,.1) 

#estimateur projection pour differents N de f3
f.func.N(data3,f3,1) 

#erreur de NW de datawave
  g = f.NW(datawave,K.gauss,hwave)
  e=datawave-g
  par(mfrow = c(1, 2))  
  plot(e,type='l')
  plot(x=g,y=e,pch=19,cex=.5)

#estimation de sigma dans une forme alternative pour data1
  data=datawave
  n = length(data)
  teta=f.estim.teta(data,n)
  s=sapply( (1:n-1), function(q){ sqrt( n * sum( (teta[(n-q):(n-1)])^2 )/q )} )
  d=density(s)
  hist(s,freq=FALSE,breaks=50,col="lightgrey",ylim=c(0,max(d$y)),main="estimation de sigma")
  points(d,type='l',lwd=2,col="red")
  