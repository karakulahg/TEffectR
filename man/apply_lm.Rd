\name{apply_lm}
\alias{apply_lm}
\title{Model gene expression levels using the given list of covariates and TE counts.}
\description{
This core function applies filtering, TMM normalization, voom transformation and LM to the given raw count expression values, respectively. It takes four arguments: (i) raw gene counts, (ii) raw TE counts, (iii) a data frame containing user-defined covariates (e.g. tissue type, disease state), and (iv) the output of get_overlaps() function. This function returns three outputs: (i) a tsv file containing the p-value of each model, significance level of covariates and associated adjusted R squared values, (ii) another tsv file containing log2(CPM) values of genes and TEs included in LM, and (iii) a group of diagnostic plots for each significant model (p < 0.05).
}
\usage{
apply_lm(gene.annotation, gene.counts, repeat.counts, covariates, prefix)
}
\arguments{
  \item{gene.annotation}{It is a data frame containing the output of get_intervals() function.}
  \item{gene.counts}{It is a data frame containing the raw read counts of genes.}
  \item{repeat.counts}{It is a data frame containing the output of summarize_repeat_counts() function.}
  \item{covariates}{It is a data frame containing the user-defined covariates (e.g. tissue type, disease state). If the TE expression as the single predictor, this parameter takes NULL value.}
  \item{prefix}{It is a string given by the user. This prefix is added to output of apply_lm() function. }

}
\seealso{
  \code{\link{}}.
}\examples{
#Create a data frame containing user-defined covariates.
df.cov<-data.frame( tissue_type=c(rep("Normal", 22), rep("Tumor", 22)), patient=c(c(1:22), c(1:22)) )

#Apply multiple linear regression models using the given list of covariates and TE counts.
TEffectR::apply_lm<-function(gene.annotation = annotation, count.matrix = count.matrix, repeat.counts = repeat.counts, covariates = df.cov, prefix="LTR-5kb-")
}

\keyword{file}
