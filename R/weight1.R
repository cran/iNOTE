weight1 <-
function(
	M, 
	G,
	C,
	Y,
	X,
	consider.gene=TRUE,
	consider.intx=TRUE
){

	## calculate MM^T, GG^T, CC^T, mu0
	MMt<-M%*%t(M)
	GGt<-G%*%t(G)
	CCt<-C%*%t(C)
	fit0<-glm(Y~X, family=binomial)
	eta0<-predict(fit0)
	mu0<-exp(eta0)/(1+exp(eta0))

	## calculate K
	kk<-mu0*(1-mu0)
	kii<--4*mu0^4+8*mu0^3-5*mu0^2+mu0
	K<-2*kk%*%t(kk)
	diag(K)<-kii	

	## calculate C3
	ci<-2*mu0^3-3*mu0^2+mu0
	C3<-diag(ci)

	## calculate W
	W<-diag(kk)

	## calculate I
	MM<-as.numeric(MMt)
	GG<-as.numeric(GGt)
	CC<-as.numeric(CCt)	
	KK<-as.numeric(K)
	
	Imm<-sum(MM*KK*MM)*0.25
	Igg<-sum(GG*KK*GG)*0.25
	Icc<-sum(CC*KK*CC)*0.25
	Img<-sum(MM*KK*GG)*0.25
	Imc<-sum(MM*KK*CC)*0.25
	Igc<-sum(GG*KK*CC)*0.25

	Iam<-t(X)%*%C3%*%diag(MMt)*0.5
	Iag<-t(X)%*%C3%*%diag(GGt)*0.5
	Iac<-t(X)%*%C3%*%diag(CCt)*0.5
	
  if(class(X)=='numeric'){
    I.aa<-as.numeric(t(X)%*%W%*%(X))
    I.tt<-matrix(c(Imm, Img, Imc, Img, Igg, Igc, Imc, Igc, Icc), 3)
    I.at<-c(Iam, Iag, Iac)
    }
	if(class(X)=='matrix'){
	  I.aa<-t(X)%*%W%*%(X)
	  I.tt<-matrix(c(Imm, Img, Imc, Img, Igg, Igc, Imc, Igc, Icc), 3)
	  I.at<-cbind(Iam, Iag, Iac)
	  }
  if(consider.gene){
  	if (!consider.intx & consider.gene & class(X)=='matrix') {I.tt<-I.tt[1:2,1:2]; I.at<-I.at[,1:2]}
    if (!consider.intx & consider.gene & class(X)=='numeric') {I.tt<-I.tt[1:2,1:2]; I.at<-I.at[1:2]}
  
  	## calculate efficient information and its inverse
    if(class(X)=='numeric'){I.til<-I.tt-I.at%*%solve(I.aa)%*%t(I.at)}
    if(class(X)=='matrix'){I.til<-I.tt-t(I.at)%*%solve(I.aa)%*%I.at}
  	a<-sqrt(1/diag(I.til))
  	a<-a/a[2]
  }
  if(!consider.gene){
    I.til<-Imm-t(Iam)%*%solve(I.aa)%*%(Iam)
    a<-sqrt(1/diag(I.til))
  }
	return(a)

}
