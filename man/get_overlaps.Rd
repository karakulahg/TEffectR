\name{get_overlaps}
\alias{get_overlaps}
\title{Getting overlapped intervals}
\description{
  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  \code{save}.
}
\usage{
  get_overlaps(g,r,strand,distance,repeat_class=NULL)
}
\arguments{
  \item{g}{}
  \item{r}{}
  \item{strand}{}
  \item{distance}{}
  \item{repeat_class}{}
}
\seealso{
  \code{\link{}}.
}\examples{
  ## getting overlapped intervals between given gene intervals and repeat intervals
  get_overlaps(g = genes, r = rm, strand = "same", distance = 10000, repeat_class = "LINE")
}

\keyword{file}
