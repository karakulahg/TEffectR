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

2. Download the most recent RepeatMasker [annotation file](http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html)

3. Read the downloaded annotation file and parse it for the downstream analysis. In our case, we use hg38 assembly:
```

repeatmasker.annotation <- TEffectR::rm_format(filepath = "~/Path2Directory/hg38.fa.out.gz" )

```
4. Read raw gene counts. An example gene count matrix can be dowloaded from: URL
```

exprs <- read.csv("gene_count_matrix.csv", row.names = 1, header=T, stringsAsFactors = F)

```
5. Retrieve the genomic locations of all genes in the given read count matrix.

    - The URL argument takes the version of Ensembl database used for gene expression quantification and can be listed [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#using-archived-versions-of-ensembl). Alternatively, you can list using the following command:  
    
        ```
        biomaRt::listEnsemblArchives()    
        ```
    
    - ID.type must be ensembl_gene_name, ensembl_gene_id, or ensembl_transcript_id.
    
    - In our case, we use Ensembl gene IDs (e.g. ENSG00000000003, ENSG00000000005, ...):
    
        ```
        gene.annotation <- get_intervals(x = rownames(exprs), assembly="hg38", ID.type = "ensembl_gene_id", URL="http://dec2014.archive.ensembl.org" ) 

        ```
6. The following function takes the genomic intervals of genes and TEs as input. Besides, the user also requires to provide three additional parameters: (i) the maximum distance allowed between the start sites of genes and TEs, (ii) whether genes and TEs must be located in same strand and (iii) TE family or subfamily name (e.g. SINE, LINE). This function helps to detect TEs that are localized within the upstream of genes of interest.

```

overlaps <- TEffectR::get_overlaps(g=gene.annotation, r=repeatmasker.annotation, strand = "same", distance = 2000, repeat_type = "LTR")

```
7. Count uniquely mapped reads to the TEs that are located within 2kb upstream of the given gene list. This step returns a raw count matrix of the total number of reads originated from TE sequences. Only the reads exhibiting 100\% overlap with given TE regions are considered and the user needs to specify individual path of each BAM file as input. All BAM files used in this step can be dowloaded from: URL This step may take up to hourse depending on the number of BAM files.
    
```

BAM.list <- c("~/Path2Directory/SRR5962198/SRR5962198_unique_sorted.bam",
             "~/Path2Directory/SRR5962199/SRR5962199_unique_sorted.bam",
             "~/Path2Directory/SRR5962200/SRR5962200_unique_sorted.bam",
             "~/Path2Directory/SRR5962201/SRR5962201_unique_sorted.bam",
             "~/Path2Directory/SRR5962202/SRR5962202_unique_sorted.bam",
             "~/Path2Directory/SRR5962203/SRR5962203_unique_sorted.bam",
             "~/Path2Directory/SRR5962204/SRR5962204_unique_sorted.bam",
             "~/Path2Directory/SRR5962205/SRR5962205_unique_sorted.bam",
             "~/Path2Directory/SRR5962206/SRR5962206_unique_sorted.bam",
             "~/Path2Directory/SRR5962207/SRR5962207_unique_sorted.bam",
             "~/Path2Directory/SRR5962208/SRR5962208_unique_sorted.bam",
             "~/Path2Directory/SRR5962209/SRR5962209_unique_sorted.bam",
             "~/Path2Directory/SRR5962210/SRR5962210_unique_sorted.bam",
             "~/Path2Directory/SRR5962211/SRR5962211_unique_sorted.bam",
             "~/Path2Directory/SRR5962212/SRR5962212_unique_sorted.bam",
             "~/Path2Directory/SRR5962213/SRR5962213_unique_sorted.bam",
             "~/Path2Directory/SRR5962214/SRR5962214_unique_sorted.bam",
             "~/Path2Directory/SRR5962215/SRR5962215_unique_sorted.bam",
             "~/Path2Directory/SRR5962216/SRR5962216_unique_sorted.bam",
             "~/Path2Directory/SRR5962217/SRR5962217_unique_sorted.bam",
             "~/Path2Directory/SRR5962218/SRR5962218_unique_sorted.bam",
             "~/Path2Directory/SRR5962219/SRR5962219_unique_sorted.bam",
             "~/Path2Directory/SRR5962220/SRR5962220_unique_sorted.bam",
             "~/Path2Directory/SRR5962221/SRR5962221_unique_sorted.bam",
             "~/Path2Directory/SRR5962222/SRR5962222_unique_sorted.bam",
             "~/Path2Directory/SRR5962223/SRR5962223_unique_sorted.bam",
             "~/Path2Directory/SRR5962224/SRR5962224_unique_sorted.bam",
             "~/Path2Directory/SRR5962225/SRR5962225_unique_sorted.bam",
             "~/Path2Directory/SRR5962226/SRR5962226_unique_sorted.bam",
             "~/Path2Directory/SRR5962227/SRR5962227_unique_sorted.bam",
             "~/Path2Directory/SRR5962228/SRR5962228_unique_sorted.bam",
             "~/Path2Directory/SRR5962229/SRR5962229_unique_sorted.bam",
             "~/Path2Directory/SRR5962230/SRR5962230_unique_sorted.bam",
             "~/Path2Directory/SRR5962231/SRR5962231_unique_sorted.bam",
             "~/Path2Directory/SRR5962232/SRR5962232_unique_sorted.bam",
             "~/Path2Directory/SRR5962233/SRR5962233_unique_sorted.bam",
             "~/Path2Directory/SRR5962234/SRR5962234_unique_sorted.bam",
             "~/Path2Directory/SRR5962235/SRR5962235_unique_sorted.bam",
             "~/Path2Directory/SRR5962236/SRR5962236_unique_sorted.bam",
             "~/Path2Directory/SRR5962237/SRR5962237_unique_sorted.bam",
             "~/Path2Directory/SRR5962238/SRR5962238_unique_sorted.bam",
             "~/Path2Directory/SRR5962239/SRR5962239_unique_sorted.bam",
             "~/Path2Directory/SRR5962240/SRR5962240_unique_sorted.bam",
             "~/Path2Directory/SRR5962241/SRR5962241_unique_sorted.bam")
             
SampleName.list<-colnames(exprs)             

TE.counts <- TEffectR::count_repeats(bamlist = BAM.list, SampleName.list, ranges=overlaps)

```

8. The following command takes the output of count_repeats() function as input. It is used to calculate the total number of sequencing reads derived from each TE that is located within the upstream of a certain gene.

```

SumOfTEs<-TEffectR::summarise_repeat_counts(counts = TE.counts, namelist = SampleName.list)

```
