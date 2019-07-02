# ###################################### <<<< 1 >>>> ###################################
# ## To Read gene or transcript expression file
x<-read.csv("~/Documents/G.Karakulah/BC/gene_count_matrix.csv", row.names = 1, header=T, stringsAsFactors = F)
x<-x[1:250,]
#
# ######################################################################################
#
#
# ###################################### <<<< 2 >>>> ###################################
#
library(Tool1)
library(stringr)
library(biomaRt)
library(biomartr)
library(dplyr)
library(Rsamtools)
#
# ## for gene annotation
#
   genes<-get_intervals(x = rownames(x), organism="hg38", ID.type = "ensembl_gene_id", URL="dec2014.archive.ensembl.org" ) #use organism parameter if you provide gene symbols
#
# ########################################################################################
#
# ###################################### <<<< 3 >>>> #####################################
# #
# # #for repeat annotation
# # download repeatmasker which you will use it from website is "http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html"
#
   rm<-rm_format(filepath="~/Downloads/hg38.fa.out.gz")
#
# #########################################################################################
#
# ###################################### <<<< 4 >>>> ######################################
# #
# #
   w<-get_overlaps(g = genes, r = rm, strand = "same", distance = 10000, repeat_family ="L1")
# #
# # write.table(w, file="overlapped.bed", quote=F, sep="\t", row.names=F, col.names=F)
# #
# #
   bamlist <- c("~/Documents/G.Karakulah/BC/uniqueBam/SRR5962198/SRR5962198_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962199/SRR5962199_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962200/SRR5962200_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962201/SRR5962201_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962202/SRR5962202_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962203/SRR5962203_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962204/SRR5962204_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962205/SRR5962205_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962206/SRR5962206_unique_sorted.bam",
                "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962207/SRR5962207_unique_sorted.bam")

   namelist<- c("SRR5962198","SRR5962199","SRR5962200","SRR5962201","SRR5962202","SRR5962203","SRR5962204","SRR5962205","SRR5962206","SRR5962207")

   e<-rm_count(bamlist = bamlist, namelist, ranges=w) #for counting

   summ<-co_summarise(e, namelist)

#
# ########################################################################################
#
#
# ######################################
#
# # for result matrix
#
#
# # df1<-genes[,5:6]
# # x<-read.csv("~/Downloads/gene_count_matrix.csv", row.names = 1, header=T, stringsAsFactors = F)
# # y<- data.frame(geneID = row.names(x), x)
# # df<-merge(df1,y,by="geneID")
# # View(df)
# # df<-df[1:4]
# # df<-df[-1]
#
#
#
# # e<-read.table("counts.txt",header = T,sep = "\t")
# # e<-e[,13:14]
# # e$ID <- seq.int(nrow(e))
# # e <- e %>% select(ID,everything())
# # colnames(e)<-colnames(df)
# # last<-rbind(df,e)
#
#
#

# # write.table(last, file="sampleMatrix1.bed", quote=F, sep="\t", row.names=F, col.names=T)
#
# # df<-read.table("sampleMatrix1.bed",header = T,sep = "\t")
# # library(edgeR)
# #
# # # DGEList object
# # norm <- DGEList(counts=df[,2:ncol(df)],genes=df[,1])
# #
# # #TMM normalization
# # norm <- calcNormFactors(norm,method = "TMM")
# #
# # normlist<-cbind(norm$genes,norm$counts)
# #
# # write.csv(normlist,"sample-TMM-Counts1.csv",row.names = F)
#
#
