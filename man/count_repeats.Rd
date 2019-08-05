\name{count_repeats}
\alias{count_repeats}
\title{Count uniquely mapped reads to the repeat regions of interest in a list of BAM files.}
\description{
This function returns a raw count matrix of the total number of reads originated from TE sequences. Only the reads exhibiting 100\% overlap with given TE regions are considered and the user needs to specify individual path of each BAM file as input.
}
\usage{
count_repeats(bamlist, namelist, ranges)
}
\arguments{
  \item{bamlist}{It is a vector which contains the paths of BAM files of interest.}
  \item{namelist}{It is a vector that contains sample names of all BAM files. Sample names must be in same order with  of them to column name to result of counts.}
  \item{ranges}{It is a Grange object. It takes the output of get_overlaps() function.}

}
\seealso{
  \code{\link{}}.
}\examples{
#Count uniquely mapped reads to the repeat regions
TEffectR::count_repeats(bamlist = bamlist, namelist = namelist, ranges = ranges)
}

\keyword{file}
