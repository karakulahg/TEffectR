\name{get_intervals}
\alias{get_intervals}
\title{Retrieve the genomic locations of all genes in a given read count matrix.}
\description{
This function is used to retrieve the genomic locations of all genes in a given read count matrix by the user. Row names of the expression matrix must be one of the following: (i) official gene symbol, (ii) Ensembl gene or (iii) transcript ID. The output of this function is utilized while determining distance between genes and TEs.
}
\usage{
get_intervals(x, assembly, ID.type, URL)
}
\arguments{
  \item{x}{a list of gene symbols, gene IDs or transcript IDs.}
  \item{assembly}{It is a character string that refers to the assembly of interest.}
  \item{ID.type}{It is a character string that refers to the ID type of interest. (ensembl_gene_name, ensemble_gene_id, ensemble_transcript_id).}
  \item{URL}{It is a character string that refers to the appropriate version of Ensembl database: http://www.ensembl.org/index.html.}
}
\seealso{
  \code{\link{}}.
}\examples{
#Retrieve genomic position of genes.
TEffectR::get_intervals(x = rownames(gene.counts), assembly="hg38", ID.type = "ensembl_gene_id", URL="dec2014.archive.ensembl.org")
}

\keyword{file}
