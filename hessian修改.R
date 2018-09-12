#Model comparison
#Part of Hessian of E_{Z|X}[F(Q,Z)]: l>k
#[this is not the upper-triangle one!]
#[reduce the for loop and exponential]
HESS<-matrix(rep(0,(K*(S+1))^2),nrow=K*(S+1))
for (j in 1:K){
  for (k in 1:S){
    #lambda
    itmp1<-numeric(N)
    itmp2<-numeric(N)
    n<-(j-1)*S+k
    a4<-sum(ALPHA[j,])
    for (i in 1:N){
      a1<-ALPHA[j,k]+X[i,k]
      a2<-sum(X[i,])+a4
      itmp1[i]<-exp(lnEZIK[i,j])*(digamma(a1)-digamma(a2)-digamma(ALPHA[j,k])+digamma(a4))
      itmp2[i]<-exp(lnEZIK[i,j])*(psigamma(a1,1)-psigamma(a2)-psigamma(ALPHA[j,k])+psigamma(a4))
    }
    HESS(n,n)<-ALPHA[j,k]*sum(itmp1)+ALPHA[j,k]*ALPHA[j,k]*sum(itmp2)
    if (k<S-1){
      for (l in k+1:S){
        n2<-(j-1)*S+l  
        itmp<-numeric(N)
        for (i in 1:N){
          a2<-sum(X[i,])+a4
          itmp[i]<-exp(lnEZIK[i,j])*(-psigamma(a2,1)+psigamma(a4,1))
        }
        HESS(n,n2)<-ALPHA[j,k]*ALPHA[j,l]*sum(itmp)
      }
    }
    #pi
    n<-K*S+k
    HESS[n,n]<- -exp(lsum(lnEZIK[,k]))/(PI[k])^2
  }
} 
