\name{rm_format}
\alias{rm_format}
\title{Formatting repeat masker file}
\description{
  Repeat annotation file is downloaded from http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html then the file is reformatted with this function.
  \code{save}.
}
\usage{
  rm_format(filepath)
}
\arguments{
  \item{filepath}{it is a character that refers directory of downloaded repeat file.}

}
\seealso{
  \code{\link{}}.
}\examples{
  ## reading and reformatting repeat masker file
   rm<-rm_format(filepath="~/Downloads/hg38.fa.out.gz")
}

\keyword{file}
