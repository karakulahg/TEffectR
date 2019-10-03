# TEffectR: An R package for studying the potential effects of transposable elements on gene expression with linear regression model
This repo is currently under review. Citation information will be provided as soon as our work is accepted. 
### What is this package used for? 
Transposable elements (TEs) are DNA sequences that are able to translocate themselves along a host genome (Biemont & Vieira 2006). This R (https://www.r-project.org) package was developed for dissecting significant associations between TEs and nearby genes in a given RNA-sequencing (RNA-seq) data set by establishing a linear regression model (LM). Our R package, namely TEffectR, makes use of publicly available RepeatMasker TE (http://www.repeatmasker.org) and Ensembl gene annotations (https://www.ensembl.org/index.html) and calculate total unique read counts of TEs from sorted and indexed genome aligned BAM files. Then, it predicts the associations of TE expressions with the transcription of adjacent genes under diverse biological conditions. These associations could be made use of by biologists to understand potential influences of TE expression on the regulation of nearby genes.

#### What are the dependencies for TEffectR ?
1. [R](https://www.r-project.org/) version should be version 3.5+
2. While using r programming, we suggest you to use [Rstudio](https://www.rstudio.com/products/rstudio/download/) which is the R statistical computing environment to use and understand functions TEffectR well.
3. [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) is required on your local computer.
4. [devtools](https://cran.r-project.org/web/packages/devtools/readme/README.html) is required to install TEffectR.
5. TEffectR uses these R packages so you have to install all of them. You may visit the following websites to install them easily: 
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

2. Download the most recent RepeatMasker [annotation file](http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html) for the organism of interest.

3. The following function takes RepeatMasker annotation file as input and extracts the genomic location of each TE along with the repeat class and family information. The output of rm_format() function is used while searching TEs overlapping in the upstream region of a given gene list. In our case, we use hg38 assembly:
```

repeatmasker.annotation <- TEffectR::rm_format(filepath = "~/Path2Directory/hg38.fa.out.gz" )

```
4. Read raw gene counts. An example gene count matrix can be dowloaded from: [here](https://drive.google.com/file/d/1icVyoqIdXqZ1jiKBAYynbTbSEl4VrtrK/view?usp=sharing). In this step, we make use of a publicly available whole transcriptome sequencing dataset including normal and tumor tissue specimens obtained from 22 ER+/HER2-breast cancer patients (GEO Accession ID: [GSE103001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103001)). The count matrix was generated with [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) - [StringTie](https://ccb.jhu.edu/software/stringtie/) pipeline. 
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
6. The following function takes the genomic intervals of genes and TEs as input. Besides, the user also requires to provide three additional parameters: (i) the maximum distance allowed between the start sites of genes and TEs, (ii) whether genes and TEs must be located in same strand and (iii) TE family or subfamily name (e.g. SINE, LINE). This function helps to detect TE species that are localized within the upstream (the "distance" argument takes positive values e.g. 5000) or downstream (the "distance" argument takes negative values e.g. -5000) of genes of interest. 

```

overlaps <- TEffectR::get_overlaps(g=gene.annotation, r=repeatmasker.annotation, strand = "strandness", distance = 5000, repeat_type = "LTR")

```
7. Count uniquely mapped reads to the TEs that are located within 5kb upstream of the given gene list. This step returns a raw count matrix of the total number of reads originated from TE sequences. Only the reads exhibiting 100\% overlap with given TE regions are considered and the user needs to specify individual paths for each BAM file as input. All BAM files used in this step can be downloaded from: [here](https://drive.google.com/file/d/1Hjac9OB07n001weLhKlYBC-GuEA5cYMv/view?usp=sharing) This step may take up to a few hours depending on the number of BAM files.
    
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
             
SampleName.list <- colnames(exprs)             

TE.counts <- TEffectR::count_repeats(bamlist = BAM.list, namelist = SampleName.list, ranges=overlaps)

```

8. The following command takes the output of count_repeats() function as input. It is used to calculate the total number of sequencing reads derived from each TE that is located within the upstream of a certain gene.

```

SumOfTEs <- TEffectR::summarize_repeat_counts(counts = TE.counts, namelist = SampleName.list)

```


9. The following core function applies filtering, TMM normalization, voom transformation and LM to the given raw count expression values, respectively. It takes four arguments: (i) raw gene counts, (ii) raw TE counts, (iii) a data frame containing user-defined covariates (e.g. tissue type, disease state), and (iv) the output of get_overlaps() function. This function returns three outputs: (i) a tsv file containing the p-value of each model, significance level of covariates and associated adjusted R squared values, (ii) another tsv file containing log2(CPM) values of genes and TEs included in LM, and (iii) a group of diagnostic plots for each significant model (p < 0.05).

```

#Create a data frame containing user-defined covariates.
df.covariates <- data.frame( tissue_type=c(rep("Normal", 22), rep("Tumor", 22)), patient=c(c(1:22), c(1:22)) )

#Apply multiple linear regression models using the given list of covariates and TE counts.
results <- TEffectR::apply_lm(gene.annotation = gene.annotation, gene.counts = exprs, repeat.counts = SumOfTEs, covariates = df.covariates, prefix = "LTR-5kb")


```

10. The tab separated file containing the p-value of each model can be downloaded [here](https://drive.google.com/file/d/1bnrhcSQTDfL2GGnUjWrImyTRcxT0_RgP/view?usp=sharing) and Log2(CPM) values of genes and repeats are available [here](https://drive.google.com/file/d/1bG_p5LagZsmohZjNWitgW0DmD97H3F6G/view?usp=sharing)


11. If you would like to run TEffectR pipeline with the TE counts generated by other quantification tools (e.g. [SQuIRE](https://www.ncbi.nlm.nih.gov/pubmed/30624635) or [TEtranscripts](https://www.ncbi.nlm.nih.gov/pubmed/26206304)), please make sure that the generated input files of such tools, which you use to run the TEffectR, must be in the following format:

```
> gene.annotation :

chr	start	end	strand	geneID	geneName
chrX	2691179	2741309		ENSG00000002586	CD99
chr4	11393150	11429765	-	ENSG00000002587	HS3ST1
chr6	136256863	136289851	-	ENSG00000029363	BCLAF1
chr14	69398015	69462388		ENSG00000029364	SLC39A9
chr3	100261000	100325251		ENSG00000036054	TBC1D23
chr1	16352575	16398145		ENSG00000055070	SZRD1
chr1	173714941	173786702		ENSG00000076321	KLHL20
chr1	45786987	46036124		ENSG00000086015	MAST2
chr18	3066807	3220108	-	ENSG00000101605	MYOM1
chr18	3247481	3256236		ENSG00000101608	MYL12A
chr19	13099033	13103161	-	ENSG00000104903	LYL1
chr11	47269161	47330031		ENSG00000110514	MADD
chr2	97655963	97664107	-	ENSG00000115073	ACTR1B
chr2	85561562	85582031		ENSG00000118640	VAMP8
chr9	100099256	100301000		ENSG00000119509	INVS
chr20	49503874	49568146	-	ENSG00000124212	PTGIS
chr20	49113339	49188367	-	ENSG00000124214	STAU1
chr20	59958427	60034011		ENSG00000124215	CDH26
chr20	50958826	50963931		ENSG00000124217	MOCS3
chrX	103062651	103064240	-	ENSG00000133169	BEX1





> gene.counts :
	SRR5962198	SRR5962199	SRR5962200	SRR5962201	SRR5962202	SRR5962203	SRR5962204	SRR5962205	SRR5962206	SRR5962207	SRR5962208
ENSG00000212040	0	0	2	0	2	0	0	1	0	0	11
ENSG00000110514	491	243	745	1233	1004	592	3493	539	412	312	1403
ENSG00000086015	438	176	702	1266	1130	526	2304	579	480	331	2122
ENSG00000211769	0	0	0	0	1	0	2	0	0	0	0
ENSG00000211768	0	0	0	0	0	0	0	0	0	0	0
ENSG00000250337	0	0	0	9	50	0	0	4	11	9	0
ENSG00000204228	31	25	69	206	125	55	760	83	59	47	81
ENSG00000255071	216	40	0	0	0	0	0	0	0	141	0
ENSG00000211767	0	0	2	0	0	0	5	0	0	0	0
ENSG00000211766	0	0	0	0	0	0	1	1	0	0	0
ENSG00000211765	0	0	0	0	0	0	0	1	0	0	0
ENSG00000211764	0	0	0	0	2	0	1	0	0	0	0
ENSG00000169740	198	166	531	723	385	365	1684	423	119	191	293
ENSG00000215869	0	0	0	0	0	0	0	0	0	0	0
ENSG00000261609	591	414	1110	2264	2144	1071	2329	568	487	449	1348
ENSG00000169744	956	497	788	1913	1344	1113	2016	566	515	534	609
ENSG00000261607	0	0	0	0	2	0	0	0	0	0	0
ENSG00000261606	0	0	0	0	0	0	0	0	6	0	22


> sum.repeat.counts :

 geneName	repeatClass	repeatName	SRR5962198	SRR5962199	SRR5962200	SRR5962201	SRR5962202	SRR5962203	SRR5962204	SRR5962205	SRR5962206	SRR5962207	SRR5962208
RPL10P1	LINE	CR1-11_Crp	0	0	0	0	0	0	0	0	0	0	0
FAM171A2	LINE	CR1-3_Croc	0	0	0	0	0	0	0	2	0	0	0
RP1-197O17.3	LINE	CR1-3_Croc	0	0	0	0	0	0	0	0	0	0	0
RP11-396O20.1	LINE	CR1-3_Croc	0	0	0	0	0	0	0	0	0	0	0
RPL7L1P5	LINE	CR1-3_Croc	0	0	0	0	0	0	0	2	0	0	0
AC108713.1	LINE	CR1-L3A_Croc	0	0	0	0	0	0	0	0	0	0	0
AC112778.1	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
ARPC2	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
CTD-2154H6.1	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
DDX10P1	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
EIF4ENIF1	LINE	HAL1	12	0	4	0	4	0	0	4	8	0	16
MIOS	LINE	HAL1	0	0	4	2	2	0	0	2	0	0	4
MIR6891	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
MOBP	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
MRPL49P1	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
MYL8P	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
PPIAP25	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
RN7SL96P	LINE	HAL1	24	0	12	0	0	0	24	0	12	0	12
RP11-109P14.2	LINE	HAL1	0	0	0	0	0	0	0	0	0	0	0
RP11-50B3.4	LINE	HAL1	0	0	8	2	26	0	4	8	0	0	24




covariates <- data.frame("TissueType" = c(rep("N",5), rep("T",6)) )

prefix<-"Test"

TEffectR::apply_lm(gene.annotation, gene.counts, repeat.counts, covariates, prefix)

```

#### Session Info

```
> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.4 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=tr_TR.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=tr_TR.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=tr_TR.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=tr_TR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] bindrcpp_0.2.2       rlist_0.4.6.1        edgeR_3.24.3         limma_3.38.3         Rsamtools_1.34.0     Biostrings_2.50.2   
 [7] XVector_0.22.0       GenomicRanges_1.34.0 GenomeInfoDb_1.18.1  IRanges_2.16.0       S4Vectors_0.20.1     BiocGenerics_0.28.0 
[13] dplyr_0.7.8          biomartr_0.8.0       biomaRt_2.38.0       stringr_1.3.1        TEffectR_0.1.0      

loaded via a namespace (and not attached):
 [1] progress_1.2.0         tidyselect_0.2.5       locfit_1.5-9.1         purrr_0.3.0            lattice_0.20-35       
 [6] yaml_2.2.0             blob_1.1.1             XML_3.98-1.16          rlang_0.3.1            pillar_1.3.1          
[11] glue_1.3.0             DBI_1.0.0              BiocParallel_1.16.5    bit64_0.9-7            GenomeInfoDbData_1.2.0
[16] bindr_0.1.1            zlibbioc_1.28.0        memoise_1.1.0          Biobase_2.42.0         curl_3.3              
[21] AnnotationDbi_1.44.0   Rcpp_1.0.0             readr_1.3.1            bit_1.1-14             hms_0.4.2             
[26] digest_0.6.18          stringi_1.2.4          grid_3.5.1             tools_3.5.1            bitops_1.0-6          
[31] magrittr_1.5           RCurl_1.95-4.11        RSQLite_2.1.1          tibble_2.0.1           crayon_1.3.4          
[36] pkgconfig_2.0.2        data.table_1.12.0      prettyunits_1.0.2      assertthat_0.2.0       httr_1.4.0            
[41] rstudioapi_0.9.0       R6_2.3.0               compiler_3.5.1        


```
