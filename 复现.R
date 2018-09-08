#[]'s need to be updated
X<-as(t(otu_table(xphy)),"matrix")
N<-nrow(X)
S<-ncol(X)
ITA<-0.1
NU<-0.1
ITER_MAX<-1000
K<-3
#EM
#initialize alpha and pi at random [see "soft k means initialiser" in source code.c?]
ALPHA<-matrix(rep(0.1,K*S),nrow = K)
PI<-runif(K,10^(-10),1-10^(-10))
EZ<-matrix(rep(0,N*K),nrow = N)
# function ln(B)
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
  ezik<-log(PI[k])+lB(ALPHA[k,]+X[i,])-lB(ALPHA[k,])-lsum(LSUMX)
  return(ezik)
}
lnEZIK<-matrix(rep(0,N*K),nrow = N)

#function for Equation 11 [to be modified!!!]
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
# Eq11<-function(lnEZIK,PI,ALPHA,X,NU,ITA){
#   EZIK<-exp(lnEZIK)
#   ezcols<-apply(EZIK,2,function(y){sum(y)})
#   lnPI<-log(PI)
#   lnBALPHA<-apply(ALPHA,1,lB)
#   E13<-sum((lnPI-lnBALPHA)*ezcols)
#   Aarr<-array(rep(ALPHA,N),dim=c(dim(ALPHA),N))
#   Xarr0<-array(rep(X,K),dim=c(dim(X),K))
#   Xarr<-aperm(Xarr0,c(3,2,1))
#   arr<-Aarr+Xarr
#   lBki<-apply(arr,c(1,3),lB)
#   E2<-sum(EZIK*t(lBki))
#   return(E13+E2-NU*sum(ALPHA)+ITA*(sum(log(ALPHA))))
# }
#update lambda by minimizing Eq11 via BFGS
#[make dtf more efficient pls!]
tf<-function(x){
  A<-matrix(exp(x),nrow=K,byrow=TRUE)
  Eq11(lnEZIK,PI,A,X,NU,ITA)
}
# dtf<-function(x){
#   der<-rep(0,K*S)
#   for (j in 1:K){
#     for (k in 1:S){
#       for (i in 1:N){
#         A1<-digamma(exp(x[(j-1)*S+k])+X[i,k])*exp(x[(j-1)*S+k])
#         sum<-0
#         for (t in 1:S){
#           sum<-sum+exp(x[(j-1)*S+t])+X[i,t]
#         }
#         A2<--digamma(sum)*exp(x[(j-1)*S+k])
#         A3<--digamma(exp(x[(j-1)*S+k]))*exp(x[(j-1)*S+k])
#         for (t in 1:S){
#           sum<-sum+exp(x[(j-1)*S+t])
#         }
#         A4<-digamma(sum)*exp(x[(j-1)*S+k])
#         der[(j-1)*S+k]<-der[(j-1)*S+k]+exp(lnEZIK[i,j])*(A1+A2+A3+A4)
#       }
#       der[(j-1)*S+k]<--der[(j-1)*S+k]+NU*exp(x[(j-1)*S+k])-ITA
#     }
#   }
#   return(der)
# }
# dtf1<-function(x){
#   x<-matrix(x,nrow=K,byrow=TRUE)
#   der<-rep(0,K*S)
#   for (j in 1:K){
#     for (k in 1:S){
#       A3<--digamma(exp(x[j,k]))*exp(x[j,k])
#       for (i in 1:N){
#         A1<-digamma(exp(x[j,k])+X[i,k])*exp(x[j,k])
#         A2<--digamma(sum((x[j,])+(X[i,])))*exp(x[j,k])
#         A4<-digamma(sum(exp(x[j,])))*exp(x[j,k])
#         der[(j-1)*S+k]<-der[(j-1)*S+k]+exp(lnEZIK[i,j])*(A1+A2+A3+A4)
#       }
#       der[(j-1)*S+k]<--der[(j-1)*S+k]+NU*exp(x[j,k])-ITA
#     }
#   }
#   return(der)
# }

dtf2<-function(x){
  x<-matrix(exp(x),nrow=K,byrow=TRUE)
  EZIK<-exp(lnEZIK)
  arows<-apply(x,1,function (y) {sum(y)})
  xrows<-apply(X,1,function (y) {sum(y)})
  matarows<-matrix(rep(arows,N),nrow=N,byrow=TRUE)
  tmp<-(digamma(matarows)+matrix(rep(xrows,K),nrow=N,byrow=FALSE))*EZIK
  A2<-matrix(rep(apply(tmp,2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)*x
  A3<-matrix(rep(apply(EZIK,2,function(y) {sum(y)}),S),nrow=K,byrow=FALSE)*digamma(x)*x
  A4<-matrix(rep(apply(EZIK*matarows,2,function(y){sum(y)}),S),nrow=K,byrow=FALSE)*x
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

lnnpi<-function(lnEZIK){
  lnNPI<-rep(0,K)
  for (kk in 1:K){
    lnNPI[k]<-lsum(lnEZIK[,k])#lnNPI[k]=ln(N*pi_k)
    return(lnNPI)
  }
}
library(gsl)
EQ11<-0
flag<-0
#[try different starting points]
for (it in 1:ITER_MAX){
  EQ11prev<-EQ11
  #update lnEZIK
  for (i in 1:N){
    for (k in 1:K){
      lnEZIK[i,k]<-LEZIK(K,PI,ALPHA,i,k)
    }
  }
  #update lambda[use different starting points!!!]
  LAMBDA0<-rep(-1,K*S)
  state<-multimin(LAMBDA0,tf,dtf2,method="bfgs")
  LAMBDA<-state$x
  #update pi
  lnNPI<-lnnpi(lnEZIK)
  }
  #[avoid this exp!!!]
  EQ11<-Eq11(lnEZIK,exp(lnNPI+log(N)),exp(LAMBDA),X,NU,ITA)
  if (abs(EQ11-EQ11prev)<0.001){
    flag<-1
    print("converged")
    break
  }
}
if (flag==0){
  print("No convergence")
}

#Model comparison
#Part of Hessian of E_{Z|X}[F(Q,Z)]: l>k
#[this is not the upper-triangle one!]
#[reduce the for loop and exponential]
HESS<-matrix(rep(0,(K*(S+1))^2),nrow=K*(S+1))
for (j in 1:K){
  for (k in 1:S){
    #lambda
    itmp1<-NULL
    itmp2<-NULL
    n<-(j-1)*S+k
    a4<-sum(ALPHA[j,])
    for (i in 1:N){
      a1<-ALPHA[j,k]+X[i,k]
      a2<-sum(X[i,])+a4
      itmp1<-c(itmp1,exp(lnEZIK[i,j])*(digamma(a1)-digamma(a2)-digamma(ALPHA[j,k])+digamma(a4)))
      itmp2<-c(itmp2,exp(lnEZIK[i,j])*(psigamma(a1,1)-psigamma(a2)-psigamma(ALPHA[j,k])+psigamma(a4)))
    }
    HESS(n,n)<-ALPHA[j,k]*sum(itmp1)+ALPHA[j,k]*ALPHA[j,k]*sum(itmp2)
    if (k<S-1){
      for (l in k+1:S){
        n2<-(j-1)*S+l  
        itmp<-NULL
        for (i in 1:N){
          a2<-sum(X[i,])+a4
          itmp<-c(itmp,exp(lnEZIK[i,j])*(-psigamma(a2,1)+psigamma(a4,1)))
        }
        HESS(n,n2)<-ALPHA[j,k]*ALPHA[j,l]*sum(itmp)
      }
    }
    #pi
    n<-K*S+k
    HESS[n,n]<- -exp(lsum(lnEZIK[,k]))/(PI[k])^2
  }
} 

