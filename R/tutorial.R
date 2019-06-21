###################################### <<<< 1 >>>> ###################################
## To Read gene or transcript expression file
# x<-read.csv("~/Documents/Kaan/gene_count_matrix.csv")  #y==hg38
# x<-read.csv("~/Documents/Kaan/transcript_count_matrix.csv") #y==hg38
# x<-scan("~/R_codes/genomeArithmetic/RepeatAnalysis/Data/genes.txt", character())

######################################################################################

###################################### <<<< 2 >>>> ###################################

# library(Tool1)

## ---> z = ensembl_transcript_id,ensembl_gene_name,ensembl_transcript_id_version,ensembl_gene_id_version
## ---> y = hg19 (Grch37) , hg38 (Grch38)
## for gene annotation
# library(stringr)
# library(biomaRt)
# names(x)[1] <-"geneID"
# genes<-a(x = x$geneID, y = "hg38", z = "ensembl_gene_id_version")

########################################################################################

###################################### <<<< 3 >>>> #####################################
#
# #for repeat annotation
# library(biomartr)
# library(dplyr)
# library(stringr)
# download repeatmasker which you will use it from website is "http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html"
# you have to format the file using the function it is  called "formatting"
#
#
#########################################################################################

###################################### <<<< 4 >>>> ######################################
#
# genes$chr<-paste("chr",genes$chr,sep = "") to edit name of chr with prefix "chr"
# w<-b(g = genes, r = repeats, strand = "same", up = 1000) for overlapping
#
# write.table(w, file="overlapped.bed", quote=F, sep="\t", row.names=F, col.names=F)
#
# e<-co(bamfilepath= "~/Documents/Kaan/NG-13693_AD_lib212351_5589_7_sorted.bam", ranges=w) for counting

########################################################################################


######################################

# for result matrix
# library(dplyr)

# > df1<-genes[,5:6]
# > df<-merge(df1,x,by="geneID")
# > View(df)
# > df<-df[-1]

# e<-read.table("counts.txt",header = T,sep = "\t")
# e<-e[,15:25]
# e$ID <- seq.int(nrow(e))
# e <- e %<% select(ID,everything())
# colnames(e)<-colnames(df)
# last<-rbind(df,e)



# write.table(last, file="sampleMatrix.bed", quote=F, sep="\t", row.names=F, col.names=T)

# df<-read.table("sampleMatrix.bed",header = T,sep = "\t")
# library(edgeR)
#
# # DGEList object
# norm <- DGEList(counts=df[,2:ncol(df)],genes=df[,1])
#
# #TMM normalization
# norm <- calcNormFactors(norm,method = "TMM")
#
# normlist<-cbind(norm$genes,norm$counts)
#
# write.csv(normlist,"sample-TMM-Counts.csv",row.names = F)


