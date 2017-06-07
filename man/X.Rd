\name{X}
\alias{X}
\docType{data}
\title{Lung Cancer Patient Clinical Data}
\description{TCGA lung cancer patient clinical data including sex, race, age at diagnosis, pathological stage of tumor at diagnosis, cancer type, and smoking history.
}
\usage{data(X)}
\format{
A named R matrix with 249 subject rows and 6 clinical phenotyp columns.
}
\details{
Clinical data were obtained from lung cancer patients (adenocarcinoma or squamous cell carcinoma) in The Cancer Genome Atlas among those who had one-year survival information, as well as both CpG methylation and mRNA expression data available. Recorded variables in this test dataset include sex, race (white or non-white), age at initial clinical diagnosis, pathological stage of tumor at biopsy, cell-type (adenocarcinoma or squamous cell), and smoking history in packyears.
}
\source{
  \describe{
  \item{TCGA}{http://cancergenome.nih.gov/}
  }
}
\examples{
data(X)
str(X)
}
\keyword{datasets}
