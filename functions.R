#--------------------------------------#
#  Creation des fonctions necessaires  #
#--------------------------------------#

#Retourne la j-ieme composante de la base de Fourrier evaluee en t
phi<-function(j,t){
  phi=c()
  for(i in j){
    if(i==1){phi=c(phi,rep(1,length(t)))}
    else if(i%%2==1){phi=c(phi,(sqrt(2)*sin(pi*(i-1)*t)))}
    else{phi=c(phi,sqrt(2)*cos(pi*i*t))}
  }
  return(phi)
}

#Retourne le vecteur teta, coordonnees de la fonction dans la base de Fourrier
f.estim.teta<-function(data,N){
  n = length(data)
  teta = sapply((1:N),function(j){mean(data*phi(j,seq(0,1,1/n)[1:n+1]))})
  
  return(teta)
}

#Renvoie la fonction estimatrice
f.estim<-function(data,N){
  teta=f.estim.teta(data,N)
  return(fun = function(t){sapply( t, function(t){sum(teta*phi((1:N),t))} )} )
}

#Compare l'erreur entre une fonction choisie et son estimee par un echantillonnage
error.function<-function(fun,data,N){
  n=length(data)
  t=seq(0,1,1/n)[1:n+1]
  return(mean((fun(t)-f.estim(data,N)(t))^2))
}

#Trace l'erreur quadratique de la fonction
plot.error.function<-function(fun,data,sigma,main){
  error<-sapply((1:length(data)-1),function(N){error.function(fun,data,N)})
  plot(error,type='l',pch=19,cex=.1,xlab='N',
       ylab='erreur quadratique',main=main,lwd=2)
  N=c(length(data)^(1/3),N.mallows(data,sigma),length(data)-1,which.min(error))
  points(x=N,y=error[N],pch=19,cex=2,col=c(2,3,4,6))
}


#Trace les differentes fonctions estimees en fonction des differents N possibles
f.func.N<-function(data,f,sigma){
  n=length(data)
  x=seq(0,1,1/n)[1:n+1]
  
  N1=as.integer(n^(1/3))
  N2=n-1
  Nopt=which.min(sapply((1:n-1),function(N){error.function(f,data,N)}))
  Nmallows=which.min(sapply((2:n-1),function(N){2*N*sigma*sigma/n-sum((f.estim.teta(data,N))^2)}))
  N=c(Nopt,N1,N2,Nmallows)
  Nnames=c('Nopt','n^(1/3)','n-1','Nmal')
  par(mfrow=c(2,2))
  for(k in (1:4)){
    plot(x,f(x),xlim=c(0,1),type='l',lwd=2,lty=2,main=paste(Nnames[k],' = ',N[k]),xlab="",ylab="")
    points(x,f.estim(data,N[k])(x),lwd=2,type='l',col='red')
    #legend("topright", legend = c('f','f.estim'), col = c(1,'red'),lty=(2:1))
  }
}

#Retourne un estimateur de sigma en fonction d'un echantillon et d'une dimension de projection q
sigma.estim<-function(data,q){
  n = length(data)
  teta=f.estim.teta(data,n)
  estim = sqrt( n*sum((teta[(n-q):(n-1)])^2) / q )
  return(estim)
}

#renvoie le graphe de l'ecart type en fonction de N
plot.sigma<-function(data,main=""){
  n = length(data)
  teta=f.estim.teta(data,n)
  s=sapply( (1:n-1), function(q){ sqrt( n * sum( (teta[(n-q):(n-1)])^2 )/q )} )
  
  plot(s,type='l',ylab ='sigma',xlab='N',lwd=2,main =main)
}


#Retourne le N de Mallows estimee pour une base de donnee et une variance connue/estimee
N.mallows<-function(data,sigma){
  n = length(data)
  teta2=f.estim.teta(data,n-1)^2
  return(which.min(sapply((1:n-1),function(N){2*N*sigma*sigma/n-sum(teta2[1:N])})))
}



#Differents noyaux utilises
K.rect<-function(u){1/2*(abs(u)<=1)}
K.gauss<-function(u){1/sqrt(2*pi)*exp(-u*u/2)}
K.trian<-function(u){(1-abs(u))*(abs(u)<=1)}
K.parab<-function(u){3/4*(1-u^2)*(abs(u)<=1)}


#Estimateur de Nadaraya-Watson pour un noyau K et une fenetre h fixe
f.NW<-function(Y,K,h){
  n=length(Y)
  t=seq(0,1,1/n)[1:n+1]
  f=c()
  for(x in t){
    S=sum(K(((1:n)/n-x)/h))
    f=c(f,sum(sapply((1:n),function(i){Y[i]*K((i/n-x)/h)}))/S)
  }
  return(f)
}

#Base de la projection
W_proj<-function(n,i,x,N){
  res=c()
  for (u in x){
    if(u !=(i/n)){
      res=c(res,1/n*sin(pi*(i/n-u)*N)/sin(pi*(i/n-u)))
    }
    else{res=c(res,N/n)}
  }
  return(res)
}

#erreur quadratique de l'estimateur de Nadaraya-Watson pour K et h fixe
error.NW<-function(fun,data,K,h){
  n=length(data)
  t=seq(0,1,1/n)[1:n+1]
  return(mean((fun(t)-f.NW(data,K,h))^2))
}


#Affiche l'erreur quadratique en fonction de la fenetre et note le hmin et h plug-in
Comp.fenetre<-function(donnees,H,fun){
  n=length(donnees)
  x=seq(0,1,1/n)[1:n+1]
  
  error.gauss<-sapply(H,function(h){error.NW(fun,donnees,K.gauss,h)})
  plot(H,error.gauss,type='l',lwd=2,xlab="h",ylab="erreur quadratique",main="Choix de h")
  #choix de la fenetre
  h=dpill(x,donnees)
  hmin=H[which.min(error.gauss)]
  points(h,error.NW(fun,donnees,K.gauss,h),col=2,pch=19,cex=2)
  points(hmin,error.NW(fun,donnees,K.gauss,hmin),col=4,pch=19,cex=2)
  legend("topleft",legend=c("h plug-in","h_opt"),col=c(2,4),pch=19)
  
}

#Affiche l'estimateur de Nadaraya-Watson pour les differents h obtenus
Comp.estim<-function(fun,data,main){
  n=length(data)
  x=seq(0,1,1/n)[-(n+1)]
  H=seq(0,0.1,0.002)
  error.gauss<-sapply(H,function(h){error.NW(fun,data,K.gauss,h)})
  
  h=dpill(x,data)
  hmin=H[which.min(error.gauss)]
  
  plot(x,data,pch=19,cex=.5,col="darkgrey",
       main=main,
       xlab="x",
       ylab="f(x)")
  
  points(x,f.NW(data,K.gauss,h),lwd=2,type='l',col=2)
  points(x,f.NW(data,K.gauss,hmin),lwd=2,type='l',col=4)
  legend("topleft",col=c(2,4),lwd=2,legend=c("NW pour h plug-in",
                                             "NW pour h optimal"))
}
