itegs <-
function(
  iCPG,           # list of J n-by-p_j CpG matrices
  iGE,            # matrix of gene expression values (n-by-J)
  iY,             # dichotomous outcome vector (n-by-1)
  iX,             # covariates (n-by-r)
  imodel='mgc',   # select disease model for whole gene set: 'mgc', 'mg', 'm', 'g'
  iapprox='pert', # select method of approximation: 'pert' or 'davies'
  i.omniseed=NA,  # J-by-1 seed vector (if using 'pert')
  no.pert=1000,   # desired # perturbations per gene level test (if using 'pert')
  gsp.emp=TRUE   # if using 'pert': return empirical iTEGS p (set to TRUE) or return approximate iTEGS p (set to FALSE) 
)
{
  gs.size=dim(iGE)[2]
  ifam=1:length(iY)
  n=length(iY)
  
  #### DAVIES METHOD ####
  if(iapprox=='davies'){
    itegs.pmin=rep(NA, gs.size)
    qhat=rep(NA, gs.size)
    
    fit.0<-glm(iY~iX, family=binomial)
    eta0<-predict(fit.0)
    mu.0<-exp(eta0)/(1+exp(eta0))
    W0.5<-diag(exp(eta0*0.5)/(1+exp(eta0)))
    R.inv=matrix(0, ncol=n, nrow=n); diag(R.inv)=1; R.inv=solve(R.inv)
    svdr<-svd(R.inv); R.inv.5<-svdr$u%*%diag(sqrt(svdr$d))

    for(g in 1:gs.size){
      if(imodel=='mgc'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                        fam=ifam, X=iX, method=iapprox, n.pert=no.pert, 
                                        consider.gene=TRUE, consider.intx=TRUE)}
      if(imodel=='mg'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                       fam=ifam, X=iX, method=iapprox, n.pert=no.pert,
                                       consider.gene=TRUE, consider.intx=FALSE)}
      if(imodel=='m'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                      fam=ifam, X=iX, method=iapprox, n.pert=no.pert,
                                      consider.gene=FALSE, consider.intx=FALSE)}
      if(imodel=='g'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                      fam=ifam, X=iX, method=iapprox, n.pert=no.pert, weight='specify',
                                      a=c(0,1,0), consider.gene=TRUE, consider.intx=FALSE)}
      
      qhat[g]=fit1$Qhat
      if(g==1){a.cent=matrix(unlist(fit1$A.cent), ncol=n)}
      if(g!=1) {a.cent=a.cent+matrix(unlist(unlist(fit1$A.cent)), ncol=n)}
    }
    
    q.geneset=sum(qhat)
    DD.MOD=eigen(W0.5%*%R.inv.5%*%a.cent%*%R.inv.5%*%W0.5/n, symmetric=TRUE)$value
    itegs.test=list(p=davies(q.geneset, lambda=DD.MOD[DD.MOD>1e-6])$Qq, iapprox='davies', imodel=imodel)
    
  }
  
  #### RESAMPLING BASED PERTURBATION ####
  if(iapprox=='pert'){
    itegs.q=rep(NA, gs.size)
    itegs.q.00=matrix(NA, ncol=gs.size, nrow=no.pert)
    
    for(g in 1:gs.size){
      if(imodel=='mgc'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                        fam=ifam, X=iX, method=iapprox, n.pert=no.pert, 
                                        consider.gene=TRUE, consider.intx=TRUE)}
      if(imodel=='mg'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                       fam=ifam, X=iX, method=iapprox, n.pert=no.pert, 
                                       consider.gene=TRUE, consider.intx=FALSE)}
      if(imodel=='m'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                      fam=ifam, X=iX, method=iapprox, n.pert=no.pert, 
                                      consider.gene=FALSE, consider.intx=FALSE)}
      if(imodel=='g'){fit1<-my.TEtest(M=as.matrix(iCPG[[g]]), G=as.matrix(iGE[,g]), Y=iY, 
                                      fam=ifam, X=iX, method=iapprox, n.pert=no.pert, weight='specify',
                                      a=c(0,1,0), consider.gene=TRUE, consider.intx=FALSE)}
      itegs.q[g]=fit1$Qhat
      itegs.q.00[,g]=fit1$QQ.0
    }
    
    itegs.q.gs=sum(itegs.q)
    itegs.q.00.gs=apply(itegs.q.00, 1, sum)
    
    if(gsp.emp==FALSE){itegs.test=list(p=qp.approx(itegs.q.gs, itegs.q.00.gs), p.emp=mean(itegs.q.gs<itegs.q.00.gs), iapprox='pert', imodel=imodel, gsp.emp=gsp.emp)
    }else{itegs.test=list(p=mean(itegs.q.gs<itegs.q.00.gs), iapprox='pert', imodel=imodel, gsp.emp=gsp.emp)}
    
  }

  return(itegs.test)
}
