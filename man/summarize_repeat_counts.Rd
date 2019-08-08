\name{summarize_repeat_counts}
\alias{summarize_repeat_counts}
\title{To summarize counts of repeats as gene name and repeat family}
\description{
Takes the output of count_repeats() function as input. It is used to calculate the total number of sequencing reads derived from each TE that is located within the upstream of a certain gene.
}
\usage{
summarize_repeat_counts(counts, namelist)
}
\arguments{
\item{counts}{It is a data frame containing the output of count_repeats() function.}
\item{namelist}{It is a vector consisting of sample names associated with BAM files. The user must provide same sample names being used in the given raw count matrix of genes. The order of sample names also must be the same as the gene count matrix.}

}
\seealso{
  \code{\link{}}.
}\examples{
#Summarize repeat counts and calculate the total number of sequencing reads derived from each TE.
TEffectR::summarize_repeat_counts(counts = repeats.counts, namelist = namelist)
}

\keyword{file}
