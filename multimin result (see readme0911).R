X<-as(t(otu_table(xphy)),"matrix")
N<-nrow(X)
S<-ncol(X)
ITA<-0.1
NU<-0.1
ITER_MAX<-1000
K<-3
ALPHA<-matrix(rep(0.001,K*S),nrow = K)
#LAMBDA<-matrix(rep(-10,K*S),nrow=K)
PI<-runif(K,10^(-10),1-10^(-10))
lnEZIK<-matrix(0,nrow=N,ncol=K)

lB<-function(x){
  mbf<-sum(lgamma(x))-lgamma(sum(x))
  return(mbf)
}
#calculate log(sum(exp(x)))
lsum<-function(x){
  m<-max(x)
  ls<-m+log(sum(exp(x-m)))
  return(ls)
}
#function for E[z_{ik}]
LEZIK<-function(K,PI,ALPHA,i,k){
  LSUMX<-numeric(K)
  for (cnt in 1:K){
    LSUMX[cnt]<-log(PI[cnt])+lB(ALPHA[cnt,]+X[i,])-lB(ALPHA[cnt,])
  }
  lnezik<-log(PI[k])+lB(ALPHA[k,]+X[i,])-lB(ALPHA[k,])-lsum(LSUMX)
  return(lnezik)
}

for(i in 1:N){
  for (k in 1:K){
    lnEZIK[i,k]<-LEZIK(K,PI,ALPHA,i,k)
  }
}
########################################################
Eq11<-function(lnEZIK,PI,ALPHA,X,NU,ITA){
  p1<-0
  for (i in 1:N){
    s<-0
    for (k in 1:K){
      s<-s+exp(lnEZIK[i,k])*(log(PI[k])+lB(ALPHA[k,]+X[i,])-lB(ALPHA[k,]))
    }
    p1<-p1+s
  }
  p2<--NU*sum(ALPHA)
  p3<-ITA*sum(log(ALPHA))
  return(-p1-p2-p3)
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

# xx<-numeric(100)
# fx<-numeric(100)
# df<-numeric(100)
# for (i in 1:100){
#   xx[i]<--5+i*0.05
#   x<-c(-10,-5+i*0.05,rep(-10,K*S-2))
#   fx[i]<-tf(x)
#   df[i]<-dtf2(x)[2]
# }
# 
# d<-numeric(100)
# plot(xx,fx,type='l')
# test<-c(10,30,60,90)
# for (j in 1:length(test)){
# for (i in 1:100){
#   tx<-xx[i]
#   d[i]<-fx[test[j]]-(xx[test[j]]-tx)*df[test[j]]
# }
#   lines(xx,d,col='red')}


x0<-rep(-10,K*S)
# state<-multimin.init(x0,tf,dtf2)
# for (i in 1:100){
#   state<-multimin.iterate(state)
#   print(c(i,state$x[10]))
# }
state<-multimin(x0,tf,dtf2,method="bfgs")
