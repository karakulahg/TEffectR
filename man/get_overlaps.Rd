\name{get_overlaps}
\alias{get_overlaps}
\title{Finding overlapping intervals}
\description{
  Finds overlaps ranges between Gene alignments and intervals of Repeats.
  \code{save}.
}
\usage{
  get_overlaps(g,r,strand,distance,repeat_class)
}
\arguments{
  \item{g}{this is a data frame which contains gene annotatitons that out of get_intervals() function.}
  \item{r}{it is a data frame which contain formatted repeat annotations that out of rm_format() function.}
  \item{strand}{it is a character string that can be "same" or "opposite". it defines strand's importance while overlapping}
  \item{distance}{it is a integer that can be negative or positive. Negative digit refers to down stream; Positive digit refers to up stream.}
  \item{repeat_class}{it is a character string that have to being selecting one of the any repeat classes}
}
\seealso{
  \code{\link{}}.
}\examples{
  ## getting overlapped intervals between given gene intervals and repeat intervals
  get_overlaps(g = genes, r = rm, strand = "same", distance = 10000, repeat_class = "LINE")
}

\keyword{file}
