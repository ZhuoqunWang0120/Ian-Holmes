# Report 4 - Reproducing Ian Holmes' method
Zhuoqun Wang
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>
## Model description and notations  
- Multinomial sampling  
Let $X$ be a matrix with elements $$X_{ij}=\text{the observed abundance of taxa j in sample i}, i=1,2,...,N, j=1,2,...,S$$, where N is the total number of samples, S the total number of taxa. Denote the $i$th row of $X$ by $\bar{X_i}$.  
Assume $\bar{X_i}$~Multinomial$\bar{p_i}$, then the likelihood of observing sample i is 
$$L_i(\bar{X_i}|\bar{p_i})=(\sum_{j=1}^S X_{ij})! \prod_{j=1}^S \frac{p_{ij}^{X_{ij}}}{X_{ij} !}$$  
The total likelihood is 
$$L(X|\bar{p_1},...,\bar{p_N}=\prod_{i=1}^N Li(X_i|\bar{p_i}))$$  
  
- Dirichlet mixture priors  
Each community vector $\bar{p_i}$ is derived from derived from a single metacommunity j with probability $\pi_{j}$. Define indicator vector $\bar{z_i}$, where $z_{ij}=1$ indicates sample i is derived from metacommunity j, then $P(\bar{z_i})=\prod_{k=1}^K \pi_k^{z_{ik}}$ where K is number of metacommunities.  
Assume $\bar{p_i}$~Dirichlet($\bar{\alpha_i}$), then the complete mixture prior is $$P(\bar{p_i}|Q)=\sum_{k=1}^K Dir(\bar(p_i)|\bar{\alpha_k})\pi_k$$  
where the mixture prior hyperparameters $Q=(K,\bar{\alpha_1},...,\bar{\alpha_K},\pi_1,...,\pi_K)$.
Let $\bar{\alpha_i}$~Gamma($\eta$,$\nu$) $(i.i.d.)$. In their source code they set $\eta=\nu=0.1$.  
  
## EM algorithm to find Q  
Reparameterisation: $\alpha_{ij}=e^{\lambda_{ij}}$  

$$\log P(X,Q)\geq E_{Z|X} P(Q,Z|X)=\sum_i\sum_k E_{Z|X}[z_{ik}](\log \pi_k+\log B(\bar{\alpha_k}+\bar{X_i})-B(\bar{\alpha_k}))-\nu\sum_{j,k} \alpha_{jk}+\eta\sum_{j,k} \log(\alpha_{jk})$$
$$\frac{\partial E_{Z|X}[F(Q,Z)]}{\partial \lambda_{jk}}=\sum_i E[z_{ij}](\Phi(\alpha_{jk}+X_{ik})\alpha_{jk}-\Phi(\sum_t (\alpha_{jt}+X_{it}))\alpha_{jk}-\Phi(\alpha_{jk})\alpha_{jk}+\Phi(\sum_t \alpha_{jt})\alpha_{jk})-\nu\alpha_{jk}+\eta$$

  
## Model comparison

$$\frac{\partial^2 E_{Z|X}[F(Q,Z)]}{\partial \lambda_{jk}^2}=\alpha_{jk}\sum_iE[z_{ij}](\Phi(\alpha_{jk}+X_{ik})-\Phi(\sum_t (\alpha_{jt}+X_{it}))-\Phi(\alpha_{jk})+\Phi(\sum_t \alpha_{jt}))+\alpha_{jk}^2(\Phi_1(\alpha_{jk}+X_{ik})-\Phi_1(\sum_t (\alpha_{jt}+X_{it}))-\Phi_1(\alpha_{jk})+\Phi_1(\sum_t \alpha_{jt}))$$  
$$\frac{\partial^2 E_{Z|X}[F(Q,Z)]}{\partial \lambda_{jk} \lambda_{jl}}=\alpha_{jk}\alpha_{jl}\sum_iE[z_{ij}](-\Phi_1(\sum_t(\alpha_{jt}+X_{it}))+\Phi_1(\sum_t \alpha_{jt}))$$  
##Import data  
We import OTU table of IBD dataset from previous file.
```r
load("D:/WANZHUQU/Desktop/HMP_IBD_16S/20180815Report.RData")
X<-as(t(otu_table(xphy)),"matrix")
```  
## Initialization  
In Ian Holmes' source code, parameters of Gamma prior for $\alpha_{ij}$ are set as $\eta=0.1,\nu=0.1$
```r
N<-nrow(X)
S<-ncol(X)
ITA<-0.1
NU<-0.1
ITER_MAX<-1000
K<-3
```
