
source("R/biomart.R")
source("R/genomicRanges.R")

# the function that to filter biomart and return genome information.
get_intervals <- function(x, organism, ID.type, URL){
  if(length(x)>0 & is.character(ID.type)){
    if(ID.type=="ensembl_gene_name"){
      df<-filterGeneName(x,organism,URL)
      return(df)
    }else if(ID.type=="ensembl_transcript_id"){
      df<-filterTranscriptID(x,organism,URL)
      return(df)
    }else if(ID.type=="ensembl_transcript_id_version"){
      df<-filterTranscriptID_V(x,organism,URL)
      return(df)
    }else if(ID.type=="ensembl_gene_id_version"){
      df<-filterGeneID_V(x,organism,URL)
      return(df)
    }else if(ID.type=="ensembl_gene_id"){
      df<-filterGeneID(x,organism,URL)
      return(df)
    }else{
      return(NULL)
    }

  }else{
    print("not found input")
    return(NULL)}

}

# the function is called as b is returned overlap positions between genes and repeats
get_overlaps <- function(g,r,strand,distance,repeat_family=NULL,repeat_class=NULL,repeat_name=NULL){
  if(is.data.frame(g) & is.data.frame(r) & is.numeric(distance)){
    if(!is.null(repeat_family)){
      hit<-r$repeat_family==repeat_family
      r<-r[hit,]
    }
    if(!is.null(repeat_class)){
      hit<-r$repeat_class==repeat_class
      r<-r[hit,]
    }
    if(!is.null(repeat_name)){
      hit<-r$repeat_name==repeat_name
      r<-r[hit,]
    }
    if(distance>0){
      g<-getUpstream(g,distance,FALSE)
    }else if(distance<0){
      g<-getDownstreams(g,distance,FALSE)
    }
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


rm_format <- function(filepath){
  dt<-read_rm(filepath)
  last<-as.data.frame(str_split_fixed(dt$matching_class, "/", 2))
  dt<-data.frame("chr"=dt$qry_id, "start"=dt$qry_start, "end"=dt$qry_end, "strand"=dt$matching_repeat, "repeat_name"=dt$repeat_id, "repeat_class"=last$V1, "repeat_family"=last$V2)
  dt$strand<-as.character(dt$strand)
  dt$strand <- replace(dt$strand, dt$strand=="C", "-")
  return(dt)
}


rm_count<-function(bamfiles,ranges){
  # counts <- count.reads( bamfiles, ranges, binary = F )
  write.table(ranges, file="overlapped.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste("bedtools multicov -s -f 1 -D -bams ", bamfiles, " -bed", "overlapped.bed > counts.txt"))
  counts<-read.csv("counts.txt",sep = "\t", header = F)
  colnames(counts)<-c(colnames(as.data.frame(ranges)),"counts")
  return(counts)
}
