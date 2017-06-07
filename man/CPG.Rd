\name{CPG}
\alias{CPG}
\docType{data}
\title{SRC Signaling Pathway CpG Methylation Data in TCGA
} 
\description{TCGA lung cancer patient Illumina 450K CpG methylation data for genes in the SRC signaling pathway geneset identified Gautschi et al (2008).
}
\usage{data(CPG)}
\format{
  A list object with eight indices, containing the matrix of CpG methylation levels for sites associated with each individual gene in the SRC Signaling pathway for 249 TCGA subjects.
  \describe{
  \item{\code{ID1}}{A named matrix with 249 subject rows and 14 CpG columns.}
  \item{\code{ID2}}{A named matrix with 249 subject rows and 15 CpG columns.}
  \item{\code{ID3}}{A named matrix with 249 subject rows and 14 CpG columns.}
  \item{\code{ID4}}{A named matrix with 249 subject rows and 8 CpG columns.}
  \item{\code{SERPINE1}}{A named matrix with 249 subject rows and 14 CpG columns.}
  \item{\code{SMAD6}}{A named matrix with 249 subject rows and 36 CpG columns.}
  \item{\code{SMAD7}}{A named matrix with 249 subject rows and 24 CpG columns.}
  \item{\code{TGFB1}}{A named matrix with 249 subject rows and 21 CpG columns.}
  }
}
\details{Genomic data were extracted from lung cancer patients (adenocarcinoma or squamous cell carcinoma) in The Cancer Genome Atlas and adjusted for batch effects. The CPG dataset represents a subset of the CpG methylation data for patients with non-missing survival status at one year after initial diagnosis in genes identified by Gautschi et al (2008) in the Molecular Signatures Database (MsigDB; Subramanian et al (2005)) as lung-cancer associated members of the SRC signaling pathway.
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
data(CPG)
str(CPG)
}
\keyword{datasets}
