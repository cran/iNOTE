\name{test.pw.info}
\alias{test.pw.info}
\docType{data}
\title{
Annotation Information for Gautschi SRC CpG Sites
}
\description{
Includes the name of the pathway within MsigDB, the full list of lung-cancer associated genes originally reported by Gautschi et al (2008) in the SRC signaling pathway, and limited annotation information for the CpG sites which both mapped to the members of the gene set and were also available for testing in the integrative analyses.
}
\usage{data(test.pw.info)}
\format{An R list object with three indices: \code{pw.name} (the name of the MsigDB pathway), \code{pw.genes} (the list of genes in the pathway), and \code{annot} (limited annotation information on the CpG sites that are associated with gene members of the pathway).
}
\source{
 \describe{
  \item{TCGA}{http://cancergenome.nih.gov/}
  \item{MsigDB}{http://software.broadinstitute.org/gsea/msigdb/}
  }
}
\references{
Gautschi O, Tepper CG, Purnell PR, Izumiya Y, Evans CP, Green TP, et al. (2008) 
Regulation of Id1 expression by SRC: implications for targeting of the bone morphogenetic protein pathway in cancer. 
\emph{Cancer Research}. \bold{68}(7):2250–8. 

Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. (2005) 
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. 
\emph{Proc Natl Acad Sci USA}. \bold{102}(43):15545–50.
}
\examples{
data(test.pw.info)
ls(test.pw.info)
}
\keyword{datasets}
