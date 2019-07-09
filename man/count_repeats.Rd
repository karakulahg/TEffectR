\name{count_repeats}
\alias{count_repeats}
\title{Counting repeats the number of overlaps in each BAM file }
\description{

  \code{save}.
}
\usage{
  count_repeats(bamlist,namelist,ranges)
}
\arguments{
  \item{bamlist}{it is a vector that contains all of bam files full directories.}
  \item{namelist}{it is a vector that refers each of them to column name to result of counts.}
  \item{ranges}{it is grange object that out of get_overlaps() function and it refers overllaped ranges between gene and repeat annotation.}

}
\seealso{
  \code{\link{}}.
}\examples{
  ## counting repeats
  count_repeats(bamlist = bamlist, namelist = namelist, ranges = ranges)
}

\keyword{file}
