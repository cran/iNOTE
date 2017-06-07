qp.approx <-
function(q, q.00){
  Q.00=q.00
  EQ=mean(Q.00)
  VQ=var(Q.00)
  kappa=VQ/(2*EQ)
  nu=2*(EQ)^2/VQ
  return(p=pchisq(q/kappa, df=nu, lower.tail=F))
}
