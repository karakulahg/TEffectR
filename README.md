# TEffectR: An R package for predicting the effects of transposable elements on gene expression with linear regression model
This repo is currently under review. Citation information will be provided as soon as our work is accepted. 
### What does this package use for? 
Transposable elements (TEs) are DNA sequences that are able to translocate themselves along a host genome (Biemont & Vieira 2006). This R (https://www.r-project.org) package, using linear regression model (LM), for dissecting significant associations between TEs and nearby genes in a given RNA-sequencing (RNA-seq) data set. Our R package, namely TEffectR, makes use of publicly available RepeatMasker TE (http://www.repeatmasker.org) and Ensembl gene annotations (https://www.ensembl.org/index.html) and calculate total unique read counts of TEs from sorted and indexed genome aligned BAM files. Then, it predicts the influence of TE expression on the transcription of adjacent genes under diverse biological conditions.

#### What are the dependencies for TEffectR ?
1. [R](https://www.r-project.org/) version should be version 3.5+
2. While using r programming we suggest you to use [Rstudio](https://www.rstudio.com/products/rstudio/download/) which is the R statistical computing environment to use and understand functions TEffectR well.
3. [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) is required on your local computer.
4. [devtools](https://cran.r-project.org/web/packages/devtools/readme/README.html) is required to install TEffectR.
5. TEffectR uses these R packages so you have to install all of them. How can you install them eaisly, please visit websites:
    - [dplyr](https://dplyr.tidyverse.org/)
    - [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
    - [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
    - [biomartr](https://cran.r-project.org/web/packages/biomartr/readme/README.html)
    - [Rsamtools](https://www.bioconductor.org/packages//2.10/bioc/html/Rsamtools.html)
    - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
    - [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
    - [rlist](https://renkun-ken.github.io/rlist/)
    - [stringr](https://github.com/tidyverse/stringr)

### How to install this R package ?
```

library(devtools)

devtools::install_github("karakulahg/TEffectR")

```

### How does it work?

1. Load the library:
```

library(TEffectR)

```

2. Download the most recent RepeatMasker annotation file(http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html)

3. Read the downloaded annotation file and parse it for the downstream analysis. In our case, we use hg38 assembly:
```

rm <- TEffectR::rm_format(filepath = "~/Downloads/hg38.fa.out.gz" )

```
4. Read raw gene counts. An example gene count matrix can be dowloaded from: URL
```

x<-read.csv("gene_count_matrix.csv", row.names = 1, header=T, stringsAsFactors = F)

```
5. Retrieve the genomic locations of all genes in the given gene read count matrix.

    - The URL option which you use annotation release is a link. you can these from [this link](https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html) or you can list by this R code:  
    
        ```
        biomaRt::listEnsemblArchives()    
        ```
    
    - ID.type can be ensembl_gene_name, ensembl_transcript_id, ensembl_transcript_id_version, ensembl_gene_id_version, ensembl_gene_id.
    
    - For this example :
    
        ```
        genes <- get_intervals(x = rownames(x), assembly="hg38", ID.type = "ensembl_gene_id", URL="dec2014.archive.ensembl.org" ) 

        ```
