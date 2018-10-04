# xglom<-tax_glom(xphy,"Genus")
readpara<-function(wd){
  para<-read.table(wd,header = FALSE,sep=",",dec=".")
  readPI<-as.numeric(para[1,-1])
  readALPHA<-t(as(para[-1,-1],"matrix"))
  para<-list()
  para$readPI<-readPI
  para$readALPHA<-readALPHA
  return(para)
}
clustratio<-function(clust,K){#Don't eliminate variable K! (Note that for K=2 you can still get all 1's)
  ratio<-numeric(K)
  for (i in 1:K){
    ratio[i]<-sum(clust==i)
  }
  ratio<-ratio/length(clust)
  return(ratio)
}
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
########################################################################
library(plotly)
otu<-otu_table(xphy)
topN<-100
TopNOTUs<-names(sort(taxa_sums(xphy),TRUE)[1:topN])
x100<-prune_taxa(TopNOTUs,xphy)
# eziklist<-list()
write.csv(data.frame(otu_table(x100)),"~/Downloads/MicrobeDMMv1.0/Data/ibd_top100.csv")
ezik10<-read.table("~/Downloads/MicrobeDMMv1.0/Data/ibd_K10_z.txt",
                  header = FALSE, sep = ",", dec = ".")
clust10<-apply(ezik10,1,function(x){which(x==max(x[2:length(x)]))-1})
sample_data(x100)$clust10<-factor(clust10)
xord <- ordinate(x100, "NMDS", "bray")
p1<-plot_ordination(x100, xord, type="samples", color="clust10") + facet_wrap(~clust10)
p1<-ggplotly(p1)
p1
parawd<-"~/Downloads/MicrobeDMMv1.0/Data/ibd_K10_mixture.txt"
ALPHA10<-readpara(parawd)$readALPHA
PI10<-readpara(parawd)$readPI
clustratio(clust10,10)-PI10
######################################################
xcd<-prune_samples(sample_data(x100)$diagnosis=="CD",x100)
xuc<-prune_samples(sample_data(x100)$diagnosis=="UC",x100)
xnon<-prune_samples(sample_data(x100)$diagnosis=="nonIBD",x100)
agecd<-sample_data(xcd)$age_at_diagnosis
ageuc<-sample_data(xuc)$age_at_diagnosis
boxplot(agecd,ageuc)
qplot(ageuc, geom = "histogram",binwidth=5) + xlab("age_UC")
qplot(agecd, geom = "histogram",binwidth=5) + xlab("age_CD")
c(min(ageuc),max(ageuc))
c(min(agecd),max(agecd))
Topnon<-names(sort(taxa_sums(xnon),TRUE)[1:5])
Topuc<-names(sort(taxa_sums(xuc),TRUE)[1:5])
Topcd<-names(sort(taxa_sums(xcd),TRUE)[1:5])
x100glom<-tax_glom(x100,"Genus")
Top10<-names(sort(taxa_sums(x100),TRUE)[1:10])
x10<-prune_taxa(Top10,x100glom)
plot_bar_arranged<-function(x,topn){
  y<-tax_glom(x,"Genus")
  Top<-names(sort(taxa_sums(y),TRUE)[1:topn])
  y<-prune_taxa(Top,y)
  topgenera<-names(sort(taxa_sums(y),TRUE))[1]
  sample_data(y)$topcount<-t(otu_table(y)[topgenera,])
  xdf<-psmelt(y)
  xdf$Genus <- factor(xdf$Genus, levels = tax_table(y)[names(sort(taxa_sums(y))),"Genus"])
  xdf$Sample<-factor(xdf$Sample,levels = rownames(sort(sample_data(y)$topcount,TRUE)))
  p <- ggplot(xdf, aes(x=Sample, y=Abundance, fill=Genus, order = as.factor(Genus)))+geom_bar(stat="identity")
  p
}
###Why are there same genus in the otu table after agglomeration?
#See below
t<-tax_table(x100glom)[,"Genus"]
table(t)
#Why are there four "uncultured"???
#因为他们来自不同的family，你是zz吗哈哈哈哈哈哈哈
