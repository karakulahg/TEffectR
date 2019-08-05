\name{count_repeats}
\alias{count_repeats}
\title{Count uniquely mapped reads to the repeat regions of interest in a given list of BAM files.}
\description{
This function returns a raw count matrix of the total number of reads originated from TE sequences. Only the reads exhibiting 100\% overlap with given TE regions are considered and the user needs to specify individual path of each BAM file as input.
}
\usage{
count_repeats(bamlist, namelist, ranges)
}
\arguments{
  \item{bamlist}{It is a vector consisting of the paths of BAM files.}
  \item{namelist}{It is a vector consisting of sample names associated with BAM files. The user must provide same sample names being used in the given raw count matrix of genes. The order of sample names also must be the same as the gene count matrix.}
  \item{ranges}{It is a GRange object. It takes the output of get_overlaps() function as input.}

}
\seealso{
  \code{\link{}}.
}\examples{
#Count uniquely mapped reads to the repeat regions.
TEffectR::count_repeats(bamlist = bamlist, namelist = namelist, ranges = ranges)
}

\keyword{file}
