\name{inote}
\alias{inote}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrative Network Omnibus Total Effect Test
}
\description{Omnibus variance-component based testing procedure to test for the total effect of methylation loci and mRNA expression across a network without requiring the specification of an underlying disease risk model.
}
\usage{
inote(iCPG, iGE, iY, iX, i.omniseed=NA, no.pert = 1000, imethod = "chi")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iCPG}{A list of J CpG matrices with dimensions nx\eqn{p_{j}}.}
  \item{iGE}{An nxJ matrix of gene expression values.}
  \item{iY}{An nx1 dichotomous outcome vector.}
  \item{iX}{An nxr covariate matrix.}
  \item{i.omniseed}{A Jx1 seed vector.}
  \item{no.pert}{No. perturbations per gene level test; defaults to 1000.}
  \item{imethod}{Omnibus testing method -- 'chi' or 'uni'}
}

\value{A list with the following components:
  \item{p}{The omnibus test p-value for the joint, integrative total effect test for the gene set.}
  \item{p.00}{The null distribution of the omnibus test statistic p-value.}
  \item{method}{The omnibus method specified by the user.}
  \item{gs.uni.mod}{If applicable, the omnibus gene-set model selected by the iNOTE-uni method.}
}
\references{
Chu S.H. and Huang Y-T. (2017)
Integrative genomic analysis of biological gene sets with applications in lung cancer. (Revise and Resubmit)
}
\author{Su H. Chu & Yen-Tsung Huang}
\examples{
data(X); data(Y); data(CPG); data(GE)
\dontrun{inote(iCPG=CPG, iGE=GE, iY=Y, iX=X, no.pert=1000, imethod='chi')}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
