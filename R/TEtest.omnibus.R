TEtest.omnibus <-
function(
  MM, 		# SNP (n-by-p matrix)
  GG, 		# gene expression (n-by-1 matrix)
  YY, 		# dichotomous outomce (n-by-1 vector)
  XX=NA,  # X covariate matrix
  fam1, 	# 1:n (n=sample size)
  nn.pert=1000, 	 # no. of perturbation
  normmix=FALSE,	 # omnibus-pvalue with gaussian mixture approximation
  omniseed=NA,
  wtMethod='pert',  # method for Qweights - 'pert' or 'davies'
  ppert.approx=TRUE
  #vtype='pert'
){
  fit1<-my.TEtest(M=MM, G=GG, Y=YY, fam=fam1, X=XX, method=wtMethod, n.pert=nn.pert,
                  consider.gene=TRUE, consider.intx=TRUE, seed=omniseed, pert.app=ppert.approx)
  fit2<-my.TEtest(M=MM, G=GG, Y=YY, fam=fam1, X=XX, method=wtMethod, n.pert=nn.pert,
                  consider.gene=TRUE, consider.intx=FALSE, seed=omniseed, pert.app=ppert.approx)
  fit3<-my.TEtest(M=MM, G=GG, Y=YY, fam=fam1, X=XX, method=wtMethod, n.pert=nn.pert, 
                  consider.gene=FALSE, consider.intx=FALSE, seed=omniseed, pert.app=ppert.approx)
  fit4<-my.TEtest(M=MM, G=GG, Y=YY, fam=fam1, X=XX, method=wtMethod, n.pert=nn.pert, weight='specify',
                  a=c(0,1,0), consider.gene=TRUE, consider.intx=FALSE, seed=omniseed, pert.app=ppert.approx)
  
  if(wtMethod=='pert'){
  
    pmin.obs<-min(fit1[[1]], fit2[[1]], fit3[[1]], fit4[[1]])
    pmin.nul<-apply(cbind(fit1[[2]], fit2[[2]], fit3[[2]], fit4[[2]]), 1, min)
    if (!normmix) pval.omb<-mean(pmin.nul<pmin.obs)
      
    if (normmix){
      mix<-normalmixEM(qnorm(pmin.nul), fast=TRUE, epsilon=1e-6)
      pi<-mix$lambda
      mu<-mix$mu
      sd<-mix$sigma
      pval.omb<-pi[1]*pnorm((qnorm(pmin.obs)-mu[1])/sd[1])+pi[2]*pnorm((qnorm(pmin.obs)-mu[2])/sd[2])
    }
      
    pval<-list(pval.omb, fit1[[1]], fit2[[1]], fit3[[1]], fit4[[1]]) #observed under all four scenarios
    names(pval)<-c("omnibus", "p_mgc", "p_mg", "p_m", "p_g")
    pvar.pert<-list(var(pmin.nul))
    names(pvar.pert)="pvar_omni"
    pweights=list(var(pmin.nul), fit1[[4]], fit2[[4]], fit3[[4]], fit4[[4]])
    names(pweights)=c("var_omni", "var_mgc", "var_mg", "var_m", "var_g")
    qhats=list(fit1[[3]], fit2[[3]], fit3[[3]], fit4[[3]])
    names(qhats)=c('Qh_mgc', 'Qh_mg', 'Qh_m', 'Qh_g')
    pnulls=list(fit1[[2]],fit2[[2]],fit3[[2]], fit4[[2]])
    names(pnulls)=c("pMGC", "pMG", "pM", 'pG')
    return(list(pval=pval, pvar.pert=pvar.pert, pmin.null=pmin.nul, 
                pnulls=pnulls,
                QQ.MGC.pert=fit1$QQ.0, QQ.MG.pert=fit2$QQ.0, QQ.M.pert=fit3$QQ.0, QQ.G.pert=fit4$QQ.0,
                pweights=pweights, Qhats=qhats))
  }
  
  if(wtMethod!='pert'){
    pval<-list(fit1[1], fit2[1], fit3[1], fit4[1])
    names(pval)=c("p_mgc", "p_mg", "p_m","p_g")
    Qhat=list(fit1[2], fit2[2], fit3[2], fit4[2])
    names(Qhat)=c('qh_mgc', 'qh_mg', 'qh_m', 'qh_g')
    A.cent<-list(fit1[3], fit2[3], fit3[3], fit4[3])
    names(A.cent)=c("a_mgc", "a_mg", "a_m", "a_g")
    return(list(pval=pval, Qhats=Qhat, A.cent=A.cent))    
  }
  
}
