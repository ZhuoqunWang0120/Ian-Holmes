#this is negative
dtf2psi1<-function(x){
  x<-matrix(exp(x),nrow=K,byrow=TRUE)
  EZIK<-exp(lnEZIK)
  arows<-apply(x,1,function (y) {sum(y)})
  xrows<-apply(X,1,function (y) {sum(y)})
  matarows<-matrix(rep(arows,N),nrow=N,byrow=TRUE)
  tmp<-(psigamma(matarows+matrix(rep(xrows,K),nrow=N,byrow=FALSE),1))*EZIK
  A2<-matrix(rep(apply(tmp,2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)*x
  A3<-matrix(rep(apply(EZIK,2,function(y) {sum(y)}),S),nrow=K,byrow=FALSE)*psigamma(x,1)*x
  A4<-matrix(rep(apply(EZIK*psigamma(matarows,1),2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)*x
  Xarr<-array(rep(X,K),dim=c(dim(X),K))
  Aarr0<-array(rep(x,N),dim=c(dim(x),N))
  Zarr0<-array(rep(EZIK,S),dim=c(dim(EZIK),S))
  Aarr<-aperm(Aarr0,c(3,2,1))#只要从dim观察即可
  Zarr<-aperm(Zarr0,c(1,3,2))
  arr<-Aarr+Xarr
  arr<-psigamma(arr,1)
  arr<-Zarr*arr
  mat<-apply(arr,2:3,function(y){sum(y)})
  A1<-t(mat)*x
  Jacobian<- -A1+A2+A3-A4+NU*x-ITA
  der<-as.vector(t(Jacobian))
  return(der)
}
#this is negative
hesspi<-function(Q){
  x<-Q[(K*S+1):(K*(S+1))]
  h<-apply(exp(lnEZIK),2,function(x){sum(x)})/x^2
  return(h)
}
#this is negative
hessdiaglam<-function(Q){
  x<-Q[1:(K*S)]
  h<-dtf2(x)+exp(x)*dtf2psi1(x)
  return(h)
}
#this is negative, 3d array,至少可以run
hesslamjkjl<-function(Q) {
  x<-Q[1:(K*S)]
  x<-matrix(exp(x),nrow=K,byrow=TRUE)
  EZIK<-exp(lnEZIK)
  arows<-apply(x,1,function (y) {sum(y)})
  xrows<-apply(X,1,function (y) {sum(y)})
  matarows<-matrix(rep(arows,N),nrow=N,byrow=TRUE)
  tmp<-(psigamma(matarows+matrix(rep(xrows,K),nrow=N,byrow=FALSE),1))*EZIK
  A2<-matrix(rep(apply(tmp,2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)
  A4<-matrix(rep(apply(EZIK*psigamma(matarows,1),2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)
  #把矩阵mat面朝前方重复n次的array：array(rep(mat,n),dim=c(dim(mat),n))
  jkjl<-array(rep((A2-A4),S),dim=c(dim(A2),S))
  Ajk<-array(rep(x,S),dim=c(dim(x),S))
  Ajl<-aperm(Ajk,c(1,3,2))
  Ajkjl<-jkjl*Ajk*Ajl
  return(Ajkjl)
}

hesslam<-function(Q){
  hjkjl<-hesslamjkjl(Q)
  hdiag<-hessdiaglam(Q)
  hdiag.vec<-as.vector(matrix(c(hdiag,rep(S,length(hdiag))),ncol=length(hdiag),byrow=TRUE))
  hesslam.vec<-as.vector(hjkjl)+hdiag.vec
  hesslam<-array(hesslam.vec,dim=dim(hjkjl))
  return(hesslam)
}

#calculate det(H)
#this calculates the det of a matrix after replacing diagonal with 0s

HDET<-function(Q){
  DETPI<-prod(hesspi(Q))
  DETLAM<-prod(apply(hesslam(Q),1,function(x){det(x)}))
  DETH<-DETPI*DETLAM
  return(DETH)
}