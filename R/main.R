# a is a function that to filter biomart and return genome information. This function gets x, y, z parameters
# x is a vector which can include gene names, gene ids or transcript ids.
# y is a string information about assembly belong to x
# z is a string information that is defined ID type ( ensembl_gene_name, gId, tId)
# x<-read.csv("~/Documents/Kaan/gene_count_matrix.csv")  #y==hg19
# x<-read.csv("~/Documents/Kaan/transcript_count_matrix.csv") #y==hg19
# x<-scan("~/R_codes/genomeArithmetic/RepeatAnalysis/Data/genes.txt", character())

source("R/biomart.R")
source("R/repeatMasker.R")
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

# data$chr <- gsub("\\d","chr\\d",)
# g is any subset of genome that is returned from called a function.
# r is a repat annotaion file
# strand is same or strandness
#up is defined upstream.
# library(dplyr)
# library(GenomicRanges)
b <- function(g,r,strand,up){
  if(is.data.frame(g) & is.data.frame(r) & is.numeric(up)){
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

library(googledrive)
# dt<-read.csv(dt$local_path,sep = "\t")
# dt$chr<-gsub("chr","",dt$chr)
download <- function(assembly){
  dt<-readRepeatMasker(assembly)
  return(dt)
}

library(seqbias)
library(Rsamtools)
#bamfilepath<-"~/Documents/Kaan/NG-13693_AD_lib212351_5589_7_sorted.bam"
c <- function(bamfilepath,ranges){
  counts <- count.reads( bamfilepath, ranges, binary = F )
  return(counts)
}
