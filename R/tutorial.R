# # ###################################### <<<< 1 >>>> ###################################
# # ## To Read gene or transcript expression file
# x<-read.csv("~/Documents/G.Karakulah/BC/gene_count_matrix.csv", row.names = 1, header=T, stringsAsFactors = F)
# x<-x[1:250,]
# x.sb<-x[7500:12500,]
# #
# # ######################################################################################
# #
# #
# # ###################################### <<<< 2 >>>> ###################################
# #

# library(stringr)
# library(biomaRt)
# library(biomartr)
# library(dplyr)
# library(Rsamtools)
# library(edgeR)
# library(rlist)
# library(limma)
#
#
# #
# # ## for gene annotation
# #
#    gene.annotation<-get_intervals(x = rownames(x), assembly="hg38", ID.type = "ensembl_gene_id", URL="dec2014.archive.ensembl.org" ) #use organism parameter if you provide gene symbols
# #
# # ########################################################################################
# #
# # ###################################### <<<< 3 >>>> #####################################
# # #
# # # #for repeat annotation
# # # download repeatmasker which you will use it from website is "http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html"
# #
#    rm<-rm_format(filepath="~/Downloads/hg38.fa.out.gz")
# #
# # #########################################################################################
# #
# # ###################################### <<<< 4 >>>> ######################################
# # #
# # #
#    overlapped.results<-get_overlaps(g = gene.annotation, r = rm, strand = "same", distance = 10000, repeat_class = "LINE")
# # #
# # # write.table(w, file="overlapped.bed", quote=F, sep="\t", row.names=F, col.names=F)
# # #
# # #
#
# bamlist <- c("~/Documents/G.Karakulah/BC/uniqueBam/SRR5962198/SRR5962198_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962199/SRR5962199_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962200/SRR5962200_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962201/SRR5962201_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962202/SRR5962202_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962203/SRR5962203_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962204/SRR5962204_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962205/SRR5962205_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962206/SRR5962206_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962207/SRR5962207_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962208/SRR5962208_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962209/SRR5962209_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962210/SRR5962210_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962211/SRR5962211_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962212/SRR5962212_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962213/SRR5962213_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962214/SRR5962214_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962215/SRR5962215_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962216/SRR5962216_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962217/SRR5962217_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962218/SRR5962218_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962219/SRR5962219_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962220/SRR5962220_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962221/SRR5962221_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962222/SRR5962222_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962223/SRR5962223_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962224/SRR5962224_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962225/SRR5962225_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962226/SRR5962226_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962227/SRR5962227_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962228/SRR5962228_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962229/SRR5962229_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962230/SRR5962230_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962231/SRR5962231_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962232/SRR5962232_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962233/SRR5962233_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962234/SRR5962234_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962235/SRR5962235_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962236/SRR5962236_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962237/SRR5962237_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962238/SRR5962238_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962239/SRR5962239_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962240/SRR5962240_unique_sorted.bam",
#              "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962241/SRR5962241_unique_sorted.bam")
#
#
# namelist<- c("SRR5962198","SRR5962199","SRR5962200","SRR5962201","SRR5962202","SRR5962203","SRR5962204","SRR5962205","SRR5962206","SRR5962207",
#              "SRR5962208","SRR5962209","SRR5962210","SRR5962211","SRR5962212","SRR5962213","SRR5962214","SRR5962215","SRR5962216","SRR5962217",
#              "SRR5962218","SRR5962219","SRR5962220","SRR5962221","SRR5962222","SRR5962223","SRR5962224","SRR5962225","SRR5962226","SRR5962227",
#              "SRR5962228","SRR5962229","SRR5962230","SRR5962231","SRR5962232","SRR5962233","SRR5962234","SRR5962235","SRR5962236","SRR5962237",
#              "SRR5962238","SRR5962239","SRR5962240","SRR5962241")

#
#    bamlist <- c("~/Documents/G.Karakulah/BC/uniqueBam/SRR5962198/SRR5962198_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962199/SRR5962199_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962200/SRR5962200_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962201/SRR5962201_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962202/SRR5962202_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962203/SRR5962203_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962204/SRR5962204_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962205/SRR5962205_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962206/SRR5962206_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962207/SRR5962207_unique_sorted.bam")
#
#    namelist<- c("SRR5962198","SRR5962199","SRR5962200","SRR5962201","SRR5962202","SRR5962203","SRR5962204","SRR5962205","SRR5962206","SRR5962207")

#    repeat.counts<-count_repeats(bamlist = bamlist, namelist, ranges= overlapped.results) #for counting
#
#   sum.repeat.counts<-summarise_repeat_counts(repeat.counts, namelist) #repeat.counts

# xc <- data.frame("Age" = c(rep(13,22), rep(24,22)),"Sex" = c(rep("F",22),rep("M",22)))

# lm.list<-apply_lm(gene.annotation = gene.annotation, gene.counts = x, repeat.counts = sum.repeat.counts, covariates = xc,"line")

