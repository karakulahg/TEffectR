\name{rm_format}
\alias{rm_format}
\title{Parse a given RepeatMasker annotation file.}
\description{
This function takes RepeatMasker annotation file as input and extracts the genomic location of each TE along with repeat class and family information. The output of rm_format() function is used while searching TEs that are located in the upstream region of the genes of interest. Repeat annotation files can be downloaded from http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html
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
#Read and parse the RepeatMasker annotation file of interest.
TEffectR::rm_format(filepath="~/FilePath/hg38.fa.out.gz")
}

\keyword{file}
