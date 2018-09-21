
library(phyloseq)
library(microbiome)
library(readxl)
library(gsl)
# wd<-"D:/WANZHUQU/Desktop/HMP_IBD_16S/taxonomic_profiles.biom"
# xim<-import_biom(wd)
# otu<-otu_table(xim)
# otu<-as(otu,"matrix")
# otumat<-otu[,-which(colnames(otu)=="206718")] 
# otumat<-otumat[,-which(colnames(otumat)=="206672")] 
# write_phyloseq(tax_table(xim), type = "all", path = getwd())
# readtax<-read_excel("D:/WANZHUQU/Desktop/HMP_IBD_16S/tax_table.xlsx")
# readtax<-as(readtax,"matrix")
# tax<-readtax[,2:7]
# rownames(tax)<-readtax[,1]
# xphy<-phyloseq(otu_table(otumat,taxa_are_rows = TRUE),tax_table(tax))
# X<-as(t(otu_table(xphy)),"matrix")  
X<-t(read.csv("C:/cygwin64/home/WANZHUQU/gsl-1.15/MicrobeDMMv1.0/MicrobeDMMv1.0/Data/Genera.csv",
                   header = TRUE, sep = ",", dec = "."))
colnames(X)<-X[1,]
X<-X[-1,]
X<-apply(X,1:2,function(y){as.numeric(y)})
N<-nrow(X)
S<-ncol(X)
ITA<-0.1
NU<-0.1
ITER_MAX<-1000
K<-4
ALPHA<-matrix(rep(1,K*S),nrow = K)
set.seed(1)
PI<-runif(K,10^(-10),1-10^(-10))
PI<-PI/sum(PI) #normalizing PI
lnEZIK<-matrix(0,nrow=N,ncol=K)
lB<-function(x){
  mbf<-sum(lgamma(x))-lgamma(sum(x))
  return(mbf)
}
lsum<-function(x){
  m<-max(x)
  ls<-m+log(sum(exp(x-m)))
  return(ls)
}
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
library(gsl)
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
  A4<-matrix(rep(apply(EZIK*digamma(matarows),2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)*x
  Xarr<-array(rep(X,K),dim=c(dim(X),K))
  Aarr0<-array(rep(x,N),dim=c(dim(x),N))
  Zarr0<-array(rep(EZIK,S),dim=c(dim(EZIK),S))
  Aarr<-aperm(Aarr0,c(3,2,1))
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
#log(N*pi)
lnnpi<-function(lnEZIK){
  lnpi<-apply(lnEZIK,2,function(x){lsum(x)})
  return(lnpi)
}
EQ11<-0
flag<-0
LAMBDA0<-rep(0,K*S)
LAMBDA<-LAMBDA0
for (it in 1:ITER_MAX){
  print(it)
  ALPHA<-matrix(exp(LAMBDA),nrow=K,byrow=TRUE)
  EQ11prev<-EQ11
  #update lnEZIK
  for (i in 1:N){
    for (k in 1:K){
      lnEZIK[i,k]<-LEZIK(K,PI,ALPHA,i,k)
    }
  }
  #update lambda
  state<-multimin(LAMBDA,tf,dtf2,method="bfgs")
  LAMBDA<-state$x
  #update pi
  lnNPI<-lnnpi(lnEZIK)
  EQ11<-Eq11(lnEZIK,exp(lnNPI-log(N)),matrix(exp(LAMBDA),nrow=K,byrow=TRUE),X,NU,ITA)
  print(c(it,EQ11prev,EQ11,EQ11prev-EQ11))
  if (abs(EQ11-EQ11prev)<0.00001){
    flag<-1
    print("converged")
    break
  }
}
if (flag==0){
  print("No convergence")
}
clust<-apply(lnEZIK,1,function(x){which(x==max(x))})