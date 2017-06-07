my.TEtest <-
function(
  
  M,			# SNP (n-by-p matrix); n=sample size, p=no. of CpGs
  G,			# gene expression (n-by-1 matrix)
  Y,			# dichotomous outcome (n-by-1 vector)
  X=NA,			  # covariates
  consider.gene=FALSE,	# (F,F)->CpG-only analysis; (T,F)->CpG/expression joint analysis (main effect only)
  consider.intx=FALSE,	# (T,T)->CpG/expression joint analysis (main effect plus interaction)
  weight="lambda",	# "lambda" or "specify"
  a=c(1,1,1),		# specify weights if weight="specify"
  R.star=NA,		
  fam=NA,			# 1:n (n=sample size)
  method="pert",		# "satt", "davies", or "pert"
  n.pert=1000,		# no. of perturbation
  pert.app=TRUE,
  seed=NA
){
  
  ## calculate n and p; generate X and C
  n<-length(Y)
  p<-dim(M)[2]
  if (class(X)=='logical') X<-rep(1, n)
  C<-matrix(NA, nrow=n, ncol=p)
  if (consider.intx) for (i in 1:n) {C[i,]<-G[i,]*M[i,]}
  
  ## standardize M, G and C
  G<-1*(G-mean(G))/sd(G)
  s.sd<-apply(M, 2, sd)
  s.m<-apply(M, 2, mean)
  M<-t((t(M)-s.m)/s.sd)
  c.sd<-apply(C, 2, sd)
  c.m<-apply(C, 2, mean)
  C<-t((t(C)-c.m)/c.sd)
  
  ## Calculate weights a
  if (weight!="specify"){
    a<-weight1(M, G, C, Y, X, consider.gene, consider.intx)
  }
    
  ## Calculate Q
  
  fit0<-glm(Y~X, family=binomial)
  eta0<-predict(fit0)
  mu.0<-exp(eta0)/(1+exp(eta0))
  if (!consider.gene & !consider.intx) A.cent<-a[1]*M%*%t(M)
  if (consider.gene & !consider.intx) A.cent<-a[1]*M%*%t(M)+a[2]*G%*%t(G)
  if (consider.intx & consider.gene) A.cent<-a[1]*M%*%t(M)+a[2]*G%*%t(G)+a[3]*C%*%t(C)
  
  ## Estimate R.star
  
  p<-dim(M)[2]
  offd<-0
  ff<-fam
  mat<-matrix(0, nrow=length(ff), ncol=length(ff))
  mat2<-matrix(0, nrow=length(ff), ncol=length(ff))
  for (i in 1:length(ff)){
    for (j in 1:length(ff)){
      if (ff[i]==ff[j]) {
        mat[i,j]<-offd
        mat2[i,j]<-offd
      }
    }
  }
  diag(mat)<-1
  diag(mat2)<-1
  R.star<-mat
  R<-mat2
  
  R.inv<-solve(R)
  W0.5<-diag(exp(eta0*0.5)/(1+exp(eta0)))
  
  
  ## davies method ##
  
  if (method=="davies"){
    
    m<-length(unique(fam))
    Q.hat<-(1/m)*t(Y-mu.0)%*%R.inv%*%A.cent%*%R.inv%*%(Y-mu.0)
    
    svdr<-svd(R.inv)
    R.inv.5<-svdr$u%*%diag(sqrt(svdr$d))
    DD<-eigen(W0.5%*%R.inv.5%*%A.cent%*%R.inv.5%*%W0.5/m, symmetric=TRUE)$value

    pval<-list(davies.p=davies(Q.hat, lambda=DD[DD>1e-6])$Qq, Qhat=Q.hat, A.cent=A.cent)
  }
  
  
  ## perturbation ##
  
  if (method=="pert"){
    
    m<-length(unique(fam))
    Q.hat<-(1/m)*t(Y-mu.0)%*%R.inv%*%A.cent%*%R.inv%*%(Y-mu.0)
    
    if (!consider.gene & !consider.intx) V<-sqrt(a[1])*M
    if ( consider.gene & !consider.intx) V<-cbind(sqrt(a[1])*M, sqrt(a[2])*G)
    if ( consider.gene &  consider.intx) V<-cbind(sqrt(a[1])*M, sqrt(a[2])*G, sqrt(a[3])*C)
    U<-cbind(X, V)
    CC<-1/m*t(U)%*%W0.5%*%R.inv%*%W0.5%*%U
    if (!is.null(dim(X))) {
      q<-dim(X)[2]
    } else {
      q<-1
    }
    Cvx<-CC[(q+1):(dim(U)[2]), 1:q]
    Cxx<-CC[1:q, 1:q]
    Av<-cbind(-Cvx%*%solve(Cxx), diag(1, dim(V)[2]))
    
    ehalf<-(1/sqrt(m))*t(U)
    QQ.0<-rep(0, n.pert)
    if(!is.na(seed)){set.seed(seed)}
    N.m<-rnorm(m*n.pert)
    N.m<-matrix(N.m, ncol=m)
    
    epsilon<-ehalf%*%((Y-mu.0)*t(N.m))
    for (r in 1:dim(Av)[1]) QQ.0<-QQ.0+(Av[r,]%*%epsilon)^2
    QQ.0<-as.numeric(QQ.0)
    
    pval.qq0<-(n.pert-rank(QQ.0)+1)/n.pert
    pval.qq0[which.max(pval.qq0)]<-1-0.5/n.pert
    
    if (!pert.app) pval<-mean(QQ.0>Q.hat[1])
    if (pert.app){
      EQ.p<-mean(QQ.0)
      VQ.p<-var(QQ.0)
      kappa.p<-VQ.p/(2*EQ.p)
      nu.p<-2*(EQ.p)^2/VQ.p
      pval<-pchisq(Q.hat[1]/kappa.p, df=nu.p, lower.tail=FALSE)
    }
    
    svdr<-svd(R.inv)
    R.inv.5<-svdr$u%*%diag(sqrt(svdr$d))    
    h<-R.inv.5%*%W0.5%*%X
    H<-h%*%solve(t(h)%*%h)%*%t(h)
    A<-(diag(1, n)-H)%*%R.inv%*%W0.5%*%A.cent%*%W0.5%*%R.inv%*%(diag(1, n)-H)
    SVQ<-2*sum(diag(A%*%R.star%*%A%*%R.star))
    
    #Pval<-list(pval, pval.qq0, Q.hat=Q.hat, Q.pert.var=var(QQ.0), QQ.0=QQ.0)
    pval<-list(pval, pval.qq0, Qhat=Q.hat, Q.pert.var=var(QQ.0), QQ.0=QQ.0, lambdaWts=a)
    
  }
  
  return(pval)
  #if (method!="pert") return(pval)
  #if (method=="pert") return(pval)
}
