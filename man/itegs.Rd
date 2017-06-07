\name{itegs}
\alias{itegs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrative Total Effect of a Gene Set Test
}
\description{
A variance-component based testing procedure to test for the total effect of methylation loci and mRNA expression across a network after specifying an underlying disease risk model which applies to all genes in the gene set of interest.}
\usage{
itegs(iCPG, iGE, iY, iX, imodel = "mgc", iapprox = "pert", 
  i.omniseed = NA, no.pert = 1000, gsp.emp = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iCPG}{A list of J CpG matrices with dimensions nx\eqn{p_{j}}.}
  \item{iGE}{An nxJ matrix of gene expression values.}
  \item{iY}{An nx1 dichotomous outcome vector.}
  \item{iX}{An nxr covariate matrix.}
  \item{imodel}{The specified disease risk model for the whole gene set: 'mgc', 'mg', or 'm'.}
  \item{iapprox}{The preferred method of gene-set p-value calculation: 'pert' or 'davies'}
  \item{i.omniseed}{A Jx1 seed vector.}
  \item{no.pert}{If using 'pert', the no. perturbations per gene level test; defaults to 1000.}
  \item{gsp.emp}{If using 'pert': The method of calculating the perturbation based p-value: empirical (TRUE) or parametric (FALSE).}
}

\value{A list with the following components:
  \item{p}{The p-value for the joint, integrative total effect test for the gene set under a pre-specified disease risk model for the whole gene set.}
  \item{iapprox}{The method of gene-set p-value calculation: 'pert' or 'davies'.}
  \item{imodel}{The user-specified disease risk model for the whole gene set.}
  \item{p.emp}{If using 'pert', and if gsp.emp is set to 'FALSE', returns the empirical p-value in addition to the approximated p-value.}
  \item{gsp.emp}{A logical value to indicate the method of calculation for the perturbation based p-value: TRUE for empirical, FALSE for parametric.}
}
\references{Chu S.H. and Huang Y-T. (2017)
Integrative genomic analysis of biological gene sets with applications in lung cancer. (Revise and Resubmit)}
\author{Su H. Chu & Yen-Tsung Huang}

\examples{
data(X); data(Y); data(CPG); data(GE)
itegs(iCPG=CPG, iGE=GE, iY=Y, iX=X, imodel='mgc', iapprox='pert', gsp.emp=FALSE);
itegs(iCPG=CPG, iGE=GE, iY=Y, iX=X, imodel='mgc', iapprox='davies');

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
