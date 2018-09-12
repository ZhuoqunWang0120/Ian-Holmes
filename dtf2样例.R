IT<-seq(10,200,10)
systime<-numeric(length(IT))
for (it in IT){
N<-it
print(it)
K<-3
S<-1000
ALPHA<-matrix(log(rep(0.5,K*S)),nrow=K)
X<-matrix(runif(N*S,0,1000), nrow = N, ncol = S)
PI<-rep(0.5,K)
lnEZIK<-matrix(0,nrow=N,ncol=K)
for(i in 1:N){
  for (k in 1:K){
    lnEZIK[i,k]<-LEZIK(K,PI,ALPHA,i,k)
  }
}
tf<-function(x){
  A<-matrix(exp(x),nrow=K,byrow=TRUE)
  Eq11(lnEZIK,PI,A,X,NU,ITA)
}

dtf2<-function(x){
  x<-matrix(exp(x),nrow=K,byrow=TRUE)
  EZIK<-exp(lnEZIK)
  arows<-apply(x,1,function (y) {sum(y)})
  xrows<-apply(X,1,function (y) {sum(y)})
  matarows<-matrix(rep(arows,N),nrow=N,byrow=TRUE)
  tmp<-(digamma(matarows+matrix(rep(xrows,K),nrow=N,byrow=FALSE)))*EZIK
  A2<-matrix(rep(apply(tmp,2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)*x
  A3<-matrix(rep(apply(EZIK,2,function(y) {sum(y)}),S),nrow=K,byrow=FALSE)*digamma(x)*x
  ############################
  A4<-matrix(rep(apply(EZIK*digamma(matarows),2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)*x
  Xarr<-array(rep(X,K),dim=c(dim(X),K))
  Aarr0<-array(rep(x,N),dim=c(dim(x),N))
  Zarr0<-array(rep(EZIK,S),dim=c(dim(EZIK),S))
  Aarr<-aperm(Aarr0,c(3,2,1))#只要从dim观察即可
  Zarr<-aperm(Zarr0,c(1,3,2))
  arr<-Aarr+Xarr
  arr<-digamma(arr)
  arr<-Zarr*arr
  mat<-apply(arr,2:3,function(y){sum(y)})
  A1<-t(mat)*x
  Jacobian<- -A1+A2+A3-A4+NU*x-ITA
  der<-as.vector(t(Jacobian))
  return(der)
}

x0<-rep(-10,K*S)

systime[it]<-system.time(state<-multimin(x0,tf,dtf2,method="bfgs"))
print(systime[it])}
# state<-multimin(x0,tf,dtf2,method="bfgs")
plot(IT,systime,type='l')