
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
get_overlaps <- function(g,r,strand,distance,repeat_class=NULL){
  if(is.data.frame(g) & is.data.frame(r) & is.numeric(distance)){
    if(!is.null(repeat_class)){
      hit<-r$repeat_class==repeat_class
      r<-r[hit,]
    }else{
      print("please choose one of the repeat options")
      return(NULL)
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


rm_count<-function(bamlist,namelist,ranges){
  bamfiles<-paste(bamlist, collapse = ' ')
  bamFile <- BamFile(bamlist[1])
  if(stringr::str_detect(seqnames(seqinfo(bamFile)),"chr")==FALSE){
    data <-as.character(lapply(seqlevels(ranges), function(x){gsub("chr", " ", x)}))
    seqlevels(ranges)<-data
  }
  df<-as.data.frame(ranges)
  df<-df[c(1,13,14,6,4,5,2,3,7,10,11,12)] # to reorder columns for bed file format
  df<-as.data.frame(apply(df,2,function(x)gsub('\\s+', '',x))) # for removing whitespaces from fields.
  df$repeat_family <- sub("^$", ".", df$repeat_family)
  write.table(df, file="overlapped.bed", quote=F, sep="\t", row.names=F, col.names=F)
  system(paste("bedtools multicov -s -f 1 -D -bams ", paste(bamfiles), " -bed", " overlapped.bed > counts.txt"))
  counts<-read.csv("counts.txt",sep = "\t", header = F)
  colnames(counts)<-c(colnames(as.data.frame(df)),namelist)
  return(counts)
}

summarise_rm_count <-function(counts,namelist){
  if(!is.null(counts) & !is.null(namelist)){
    col_indexes <- which(colnames(counts) %in% namelist)
    b<-aggregate(list(counts[,col_indexes]), by=list(geneName=counts$geneName, repeatClass=counts$repeat_class ,repeatFamily=counts$repeat_family), FUN=sum)
    return(b)
  }
}

merge_counts<-function(gene.annotation, genes.expr,repeats.expr){
  df1<-genes[,5:6]
  y<- data.frame(geneID = row.names(genes.expr), genes.expr)
  df<-merge(df1,y,by="geneID")
  df<-df[,2:ncol(df)]

  e1<-repeats.expr[,4:ncol(repeats.expr)]
  e1$geneName <- seq.int(nrow(e1))
  e1 <- e1 %>% select(geneName,everything())
  last<-rbind(df,e1)
  return(last)
}

get_dge_transformation<- function(count.matrix){
  dge<-DGEList(count.matrix[,2:ncol(count.matrix)])
  keep <- filterByExpr(dge, min.total.count=10)
  dge <- dge[keep, keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge)
  return(v)
}


apply_lm<-function(dge.transformation,repeat.counts,group){
  l<-list()
  gene.counts<-count.matrix[row.names(dge$E),]
  s3<-group
  if(length(s3)==(ncol(gene.counts)-1)){
    for (r in 1:nrow(gene.counts)) {
      s1<-gene.counts[r,2:ncol(gene.counts)]
      hit<-repeat.counts$geneName==gene.counts$geneName[r]
      s2<-repeat.counts[hit,]
      if(nrow(s2)>0){
        goalsMenu <- s2$repeatFamily
        output <- as.data.frame(matrix(rep(0, 2 + length(goalsMenu)), nrow=1))
        names(output) <- c(gene.counts$geneName[r], as.vector(s2$repeatFamily),"group")
        output[1:length(s1),1]<-as.numeric(s1)
        for (var in 1:(ncol(output)-2)){
          output[1:length(s2[var,4:ncol(s2)]),var+1]<- as.numeric(s2[var,4:ncol(s2)])
        }
        output[1:length(s3),ncol(output)]<-as.character(s3)
        colnames(output)[1]<-gsub("-","_",colnames(output)[1])
        if(r>1){
          formula <- paste(as.character(colnames(output)[1]),"~.",sep = "")
          lm.out<-lm(formula = formula , data=output)
          l<-list.append(l,lm.out)
          # l<-list.append(l,output)
        }else if (r==1){
          formula <- paste(as.character(colnames(output)[1]),"~.",sep = "")
          lm.out<-lm(formula = formula , data=output)
          l<-list.append(list(lm.out))
          # l<-list.append(list(output))
        }
      }

    }
    return(l)
  }else{
    print("number of group does not match with sample number of gene count matrix please check it ! ")
    return(NULL)
  }

}




