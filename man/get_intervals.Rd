\name{get_intervals}
\alias{get_intervals}
\title{Getting gene intervals}
\description{
    Gene Annotations is filtered as a group of gene id or transcript id from Ensembl Genome Browser with using BiomaRt R package.
  \code{save}.
}
\usage{
get_intervals(x, organism, ID.type, URL)
}
\arguments{
  \item{x}{a rownames vector of gene.counts or transcript.counts matrix that it includes gene ids or transcript ids.}
  \item{organism}{it is a character string that defines assembly which working on.}
  \item{ID.type}{it is a character string that defines id type for biomaRt using (ensemble_gene_id, ensemble_transcript_id).}
  \item{URL}{the url is defined to interest on spesific ensemble archive version to filter gene interval (you can select the path from here : http://www.ensembl.org/index.html).}
}
\seealso{
  \code{\link{}}.
}\examples{
## getting gene intervals
get_intervals(x = rownames(gene.counts), organism="hg38", ID.type = "ensembl_gene_id", URL="dec2014.archive.ensembl.org")
}

\keyword{file}
