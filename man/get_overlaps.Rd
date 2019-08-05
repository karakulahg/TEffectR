\name{get_overlaps}
\alias{get_overlaps}
\title{Find overlapping intervals.}
\description{
Takes the genomic intervals of genes and TEs as input. Besides, the user also requires to provide three additional parameters: (i) the maximum distance allowed between the start sites of genes and TEs, (ii) whether genes and TEs must be located in same strand and (iii) TE family or subfamily name (e.g. SINE, LINE). This function helps to detect TEs that are localized within the upstream of genes of interest.
}
\usage{
get_overlaps(g, r, strand, distance, repeat_type)
}
\arguments{
  \item{g}{It is a data frame containing gene annotations retrieved via get_intervals() function.}
  \item{r}{It is a data frame containing reformatted repeat annotations created via rm_format() function.}
  \item{strand}{This argument takes one of the two values: either "same" or "strandness". In case "same" option is used, gene and TE must be in same strand.}
  \item{distance}{It is a integer value that can be negative or positive. This argument is used for determining the maximum distance allowed between the start sites of genes and TEs.}
  \item{repeat_type}{This argument takes a character string and it is used to filter a specific repeat type for the downstream analysis.}
}
\seealso{
  \code{\link{}}.
}\examples{
#Find TEs that are located within the upstream of the given gene list.
TEffectR::get_overlaps(g = genes, r = rm, strand = "same", distance = 1000, repeat_type = "LINE")
}

\keyword{file}
