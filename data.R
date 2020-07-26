#---------------------------------------#
#----- Creation de jeu de donnees  -----#
#---------------------------------------#

#Fonction 1
n=100
x=seq(0,1,1/n)[1:n+1]
f1<-function(x){(1-x^2)*sin(pi*x)*(x>1/2)}
data1=f1(x)+rnorm(n,0,0.1)

#Fonction 2
n=100
x=seq(0,1,1/n)[1:n+1]
f2<-function(x){sin(pi*x/3)+cos(pi*x*2)}
data2=f2(x)+rnorm(n,0,0.1)

#Fonction 3
n=250
x=seq(0,1,1/n)[1:n+1]
f3<-function(x){sapply(x,function(u){min(10,1/abs(cos(4*pi*u)))})}
data3=f3(x)+rnorm(n,0,1)

#Fonction 4
n=200
x=seq(0,1,1/n)[1:n+1]
f4<-function(x){4*sin(4*pi*x)+3*cos(6*pi*x)-2*sin(6*pi*x)}
data4=f4(x)+rnorm(n,0,3)

#donnees nottem
data=as.vector(nottem)

#Wave data
#https://lambda.gsfc.nasa.gov/product/map/dr1/map_tt_powspec.cfm
datawave <- read.table("dataCMBwave.txt", row.names=1, quote="\"")$V2[1:700]
