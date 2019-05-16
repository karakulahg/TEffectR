# a is a function that to filter biomart and return genome information. This function gets x, y, z parameters
# x is a vector which can include gene names, gene ids or transcript ids.
# y is a string information about assembly belong to x
# z is a string information that is defined ID type ( ensemble_gene_name, gId, tId)
#x<-read.csv("~/Documents/Kaan/gene_count_matrix.csv")  #y==hg19
#x<-read.csv("~/Documents/Kaan/transcript_count_matrix.csv") #y==hg19
#x<-scan("~/R_codes/genomeArithmetic/RepeatAnalysis/Data/genes.txt", character())
source("R/biomart.R")
source("R/repeatMasker.R")
source("R/genomicRanges.R")

a<-function(x,y,z){

if(!is.na(x) & !is.na(z)){
  if(z=="ensemble_gene_name"){
    df<-filterGeneName(x,y)
    return(df)
  }else if(z=="ensembl_transcript_id"){
    df<-filterTranscriptID(x,y)
    return(df)
  }else if(z=="ensembl_gene_id"){
    df<-filterGeneID(x,y)
    return(df)
  }

}else{
  print("not found input")
  return(NULL)}

}

# g is any subset of genome that is returned from called a function.
# assembly
# strand is same or strandness
#up is defined upstream.

b<-function(g,assembly,strand,up){
  r<-readRepeatMasker(assembly)
  # if(!is.na(g) & !is.na(r)){
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
  # }else{
  #   print("not found input")
  #   return(NULL)}
}

