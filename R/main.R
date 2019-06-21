# a is a function that to filter biomart and return genome information. This function gets x, y, z parameters
# x is a vector which can include gene names, gene ids or transcript ids.
# y is a string information about assembly belong to x
# z is a string information that is defined ID type ( ensembl_gene_name, gId, tId)
# x<-read.csv("~/Documents/Kaan/gene_count_matrix.csv")  #y==hg19
# x<-read.csv("~/Documents/Kaan/transcript_count_matrix.csv") #y==hg19
# x<-scan("~/R_codes/genomeArithmetic/RepeatAnalysis/Data/genes.txt", character())

source("R/biomart.R")
source("R/genomicRanges.R")

a <- function(x,y,z){
  if(length(x)>0 & is.character(z)){
    if(z=="ensembl_gene_name"){
      df<-filterGeneName(x,y)
      return(df)
    }else if(z=="ensembl_transcript_id"){
      df<-filterTranscriptID(x,y)
      return(df)
    }else if(z=="ensembl_transcript_id_version"){
      df<-filterTranscriptID_V(x,y)
      return(df)
    }else if(z=="ensembl_gene_id_version"){
      df<-filterGeneID_V(x,y)
      return(df)
    }else if(z=="ensembl_gene_id"){
      df<-filterGeneID(x,y)
      return(df)
    }else{
      return(NULL)
    }

  }else{
    print("not found input")
    return(NULL)}

}

# the function is called as b is returned overlap positions between genes and repeats
# data$chr <- gsub("\\d","chr\\d",)
# g is any subset of genome that is returned from called a function.
# r is a repat annotaion file
# strand is same or strandness
# up is defined upstream.
# library(dplyr)
# library(GenomicRanges)
b <- function(g,r,strand,up,family,class){
  if(is.data.frame(g) & is.data.frame(r) & is.numeric(up)){
    if(!is.na(family)){
      hit<-r$repeat_family==family
      r<-r[hit,]
    }
    if(!is.na(class)){
      hit<-r$repeat_class==class
      r<-r[hit,]
    }
    g<-getUpstream(g,up,FALSE)
    if(strand=="same"){
      g<-makeGRangeObj(g)
      r<-makeGRangeObj(r)
      overlaps<-toFindOverlaps(r,g)
      return(overlaps)
    }else if(strand=="strandness"){
      g<-makeGrObj_Unstrand(g)
      r<-makeGrObj_Unstrand(r)
      overlaps<-toFindOverlaps(r,g)
      return(overlaps)
    }
  }else{
    print("not found input")
    return(NULL)}
}


formatting <- function(filepath){
  dt<-read_rm(filepath)
  last<-as.data.frame(str_split_fixed(dt$matching_class, "/", 2))
  dt<-data.frame("chr"=dt$qry_id, "start"=dt$qry_start, "end"=dt$qry_end, "strand"=dt$matching_repeat, "repeat"=dt$repeat_id, "repeat_class"=last$V1, "repeat_family"=last$V2)
  dt$strand<-as.character(dt$strand)
  dt$strand <- replace(dt$strand, dt$strand=="C", "-")
  return(dt)
}

#bamfiles<-"~/Documents/Kaan/NG-13693_AD_lib212351_5589_7_sorted.bam"
co<-function(bamfiles,ranges){
  # counts <- count.reads( bamfilepath, ranges, binary = F )
  write.table(ranges, file="overlapped.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste("bedtools multicov -bams ", bamfiles, " -bed", "overlapped.bed > counts.txt"))
  counts<-read.csv("counts.txt",sep = "\t", header = F)
  colnames(counts)<-c(colnames(as.data.frame(ranges)),"counts")
  return(counts)
}
