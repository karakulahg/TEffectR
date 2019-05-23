# x<-read.csv("~/Documents/Kaan/gene_count_matrix.csv")  #y==hg19
# x<-read.csv("~/Documents/Kaan/transcript_count_matrix.csv") #y==hg19
# x<-scan("~/R_codes/genomeArithmetic/RepeatAnalysis/Data/genes.txt", character())
#
#
# library(Tool1)
# #ensembl_transcript_id,ensembl_gene_name,ensembl_transcript_id_version,ensembl_gene_name_version
# #for gene annotation
# genes<-a(x = x, y = "hg19", z = "ensembl_gene_name")
#
#
#
# #for repeat annotation
# library(googledrive)
# dt<-download("hg19")
# repeats<-read.csv(dt$local_path,sep = "\t")
#
# w<-b(g = genes, r = repeats, strand = "same", up = 1000)

# e<-co(bamfilepath= "~/Documents/Kaan/NG-13693_AD_lib212351_5589_7_sorted.bam", ranges=w)

