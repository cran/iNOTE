inote <-
function(
  iCPG,           # list of J n-by-p_j CpG matrices
  iGE,            # matrix of gene expression values (n-by-J)
  iY,             # dichotomous outcome vector (n-by-1)
  iX,             # covariates (n-by-r)
  i.omniseed=NA, # J-by-1 seed vector
  no.pert=1000,  # desired # perturbations per gene level test
  imethod='chi'  # 'chi' or 'uni' 
)
  {
  gs.size=dim(iGE)[2]
  if(is.na(i.omniseed)){i.omniseed=sample(1:(1000*gs.size), gs.size)}
  fam=1:length(iY)
  
  ## chi ##
  if(imethod=='chi'){
    chi.pmin.omni=rep(NA, gs.size)
    chi.pmin.omni.00=matrix(NA, ncol=gs.size, nrow=no.pert)
    
    #calculate gene level joint statistics
    for(g in 1:gs.size){
      testGene=TEtest.omnibus(MM=as.matrix(iCPG[[g]]), GG=as.matrix(iGE[,g]), 
                              YY=iY, XX=iX, fam1=1:length(iY), omniseed=i.omniseed[g], 
                              nn.pert=no.pert, ppert.approx=TRUE)
      
      chi.omni.g.p=min(unlist(testGene$pval)[-1])
      
      #approximate p-values if chi.omni.g.p==0 for gene level p value
      if(chi.omni.g.p==0){
        index.zero=which.min(testGene$pval[-1])
        min.Q.hat=unlist(testGene$Qhat)[index.zero]
        min.Q.00=unlist(list(testGene$QQ.SGC.pert, testGene$QQ.SG.pert, testGene$QQ.S.pert)[index.zero])
        chi.omni.g.p=qp.approx(t$min.Q.hat, t$min.Q.00)}
      
      # calculate gene level chi-transformed min p statistic for each gene, and its p-value
      chi.pmin.omni[g]=qchisq(chi.omni.g.p, 1, ncp = 0, lower.tail=F)
      chi.pmin.omni.00[,g]=qchisq(unlist(testGene$pmin.null), df=1, lower.tail=F)
    }
    
    # obtain network level test statistic and empirical distribution
    inote.test=list(p=mean(sum(chi.pmin.omni)<apply(chi.pmin.omni.00,1,sum)),
                    p.00=apply(chi.pmin.omni.00, 1, sum), method='chi')
  }
  
  ## uni ##
  if(imethod=='uni'){
    qhat=matrix(NA, nrow=gs.size, ncol=4)
    q00.mgc=matrix(NA, ncol=gs.size, nrow=no.pert)
    q00.mg=matrix(NA, ncol=gs.size, nrow=no.pert)
    q00.m=matrix(NA, ncol=gs.size, nrow=no.pert)
    q00.g=matrix(NA, ncol=gs.size, nrow=no.pert)
    
    #obtain gene level joint statistics
    for(g in 1:gs.size){
      testGene=TEtest.omnibus(MM=as.matrix(iCPG[[g]]), GG=as.matrix(iGE[,g]), 
                              YY=iY, XX=iX, fam1=1:length(iY), omniseed=i.omniseed[g], 
                              nn.pert=no.pert, ppert.approx=TRUE)
      
      q00.mgc[,g]=unlist(testGene$QQ.MGC.pert)
      q00.mg[,g]=unlist(testGene$QQ.MG.pert)
      q00.m[,g]=unlist(testGene$QQ.M.pert)
      q00.g[,g]=unlist(testGene$QQ.G.pert)
      
      qhat[g,]=unlist(testGene$Qhats)
    }
    
    # combine gene level joint statistics
    q.geneset=apply(qhat,2,sum)
    q00.mgc.gs=apply(q00.mgc,1,sum)
    q00.mg.gs=apply(q00.mg,1,sum)
    q00.m.gs=apply(q00.m,1,sum)
    q00.g.gs=apply(q00.g,1,sum)
    
    # calculate p-values of empirical null distributions for each model
    q00.mgc.gs.p=(no.pert-rank(q00.mgc.gs)+1)/no.pert
    q00.mg.gs.p=(no.pert-rank(q00.mg.gs)+1)/no.pert
    q00.m.gs.p=(no.pert-rank(q00.m.gs)+1)/no.pert
    q00.g.gs.p=(no.pert-rank(q00.g.gs)+1)/no.pert
    
    # select omnibus statistic and obtain p-value
    p00.gs=apply(cbind(q00.mgc.gs.p, q00.mg.gs.p, q00.m.gs.p,q00.g.gs.p), 1, min)
    gs.p=c(MGC=mean(q.geneset[1]<q00.mgc.gs), MG=mean(q.geneset[2]<q00.mg.gs), M=mean(q.geneset[3]<q00.m.gs), G=mean(q.geneset[4]<q00.g.gs))
    
    inote.test=list(p=mean(min(gs.p)>p00.gs),
                    p.00=p00.gs, gs.uni.mod=names(which.min(gs.p)), method='uni')
  }
  
  return(inote.test)
}
