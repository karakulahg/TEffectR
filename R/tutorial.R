# x<-read.csv("~/Documents/Kaan/gene_count_matrix.csv")  #y==hg19
# x<-read.csv("~/Documents/Kaan/transcript_count_matrix.csv") #y==hg19
# x<-scan("~/R_codes/genomeArithmetic/RepeatAnalysis/Data/genes.txt", character())


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



# > names(x)[1] <-"geneID"
# > df1<-genes[,5:6]
# > df<-merge(df1,x,by="geneID")
# > View(df)
# > df<-df[-1]

# genes$chr<-paste("chr",genes$chr,sep = "")

# write.table(df, file="sampleMatrix.bed", quote=F, sep="\t", row.names=F, col.names=T)

# df<-read.table("sampleMatrix.bed",header = T,sep = ",")
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


