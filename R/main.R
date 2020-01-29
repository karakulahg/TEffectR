
source("R/biomart.R")
source("R/genomicRanges.R")


# the function that to filter biomart and return genome information.
get_intervals <- function(x, assembly, ID.type, URL){
  if(length(x)>0 & is.character(ID.type)){
    if(ID.type=="ensembl_gene_name"){
      df<-filterGeneName(x,assembly,URL)
      return(df)
    }else if(ID.type=="ensembl_transcript_id"){
      df<-filterTranscriptID(x,assembly,URL)
      return(df)
    }else if(ID.type=="ensembl_transcript_id_version"){
      df<-filterTranscriptID_V(x,assembly,URL)
      return(df)
    }else if(ID.type=="ensembl_gene_id_version"){
      df<-filterGeneID_V(x,assembly,URL)
      return(df)
    }else if(ID.type=="ensembl_gene_id"){
      df<-filterGeneID(x,assembly,URL)
      return(df)
    }else{
      return(NULL)
    }

  }else{
    print("not found input")
    return(NULL)}

}

# the function is called as b is returned overlap positions between genes and repeats
get_overlaps <- function(g,r,strand,distance,repeat_type){
  options(warn=-1)
  if(is.data.frame(g) & is.data.frame(r) & is.numeric(distance)){
    if(!is.null(repeat_type)){
      hit<-r$repeat_type==repeat_type
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

    g<-makeGRangeObj(g)
    r<-makeGRangeObj(r)
    overlaps<-toFindOverlaps(r,g,strand)
    return(overlaps)

  }else{
    print("not found input")
    return(NULL)}
}


rm_format <- function(filepath){
  dt<-biomartr::read_rm(filepath)
  last<-as.data.frame(stringr::str_split_fixed(dt$matching_class, "/", 2))
  dt<-data.frame("chr"=dt$qry_id, "start"=dt$qry_start, "end"=dt$qry_end, "strand"=dt$matching_repeat, "repeat_name"=dt$repeat_id, "repeat_type"=last$V1, "repeat_family"=last$V2)
  dt$strand<-as.character(dt$strand)
  dt$strand <- replace(dt$strand, dt$strand=="C", "-")
  return(dt)
}


count_repeats<-function(bamlist,namelist,ranges){
  bamfiles<-paste(bamlist, collapse = ' ')
  bamFile <- Rsamtools::BamFile(bamlist[1])
  if(stringr::str_detect(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(bamFile)),"chr")==FALSE){
    data <-as.character(lapply(GenomeInfoDb::seqlevels(ranges), function(x){gsub("chr", " ", x)}))
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

summarize_repeat_counts <-function(counts,namelist){
  if(!is.null(counts) & !is.null(namelist)){
    col_indexes <- which(colnames(counts) %in% namelist)
    b<-aggregate(list(counts[,col_indexes]), by=list(geneName=counts$geneName, repeatClass=counts$repeat_type ,repeatName=counts$repeat_name), FUN=sum)
    return(b)
  }
}


apply_lm<-function(gene.annotation, gene.counts, repeat.counts, covariates=NULL, prefix){

  # to merge counts
  df1<-gene.annotation[,5:6]
  y<- data.frame(geneID = row.names(gene.counts), gene.counts)
  df<-merge(df1,y,by="geneID") # finding gene names as geneid
  df<-df[,2:ncol(df)]
  e1<-repeat.counts[,4:ncol(repeat.counts)]
  e1$geneName <- seq.int(nrow(e1))
  e1 <- e1 %>% select(geneName,everything())
  count.matrix<-rbind(df,e1)


  # to apply filter and TMM normalization then apply voom
  dge <- edgeR::DGEList(count.matrix[,2:ncol(count.matrix)])
  keep <- edgeR::filterByExpr(dge, min.total.count=10)
  dge <- dge[keep, keep.lib.sizes=FALSE]
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge)


  vall<-count.matrix[row.names(v$E),]
  new_voom<-cbind(data.frame("geneName"=c(vall$geneName)),v$E)
  v_ids_for_repeats<-setdiff(vall[,1],df$geneName) #to get row ids of repeats which are passed from voom translation
  col_indexes <- which(vall[,1] %in% v_ids_for_repeats)
  temp<-repeat.counts[v_ids_for_repeats,] # repeat counts with assoiated with genes
  tt<-paste(temp$geneName,temp$repeatClass,temp$repeatName,sep = ":")
  vall$geneName[col_indexes]<-tt


  prepLMdata_repeats<-temp
  prepLMdata_repeats[,4:ncol(prepLMdata_repeats)]<-new_voom[col_indexes,2:ncol(new_voom)]

  v_ids_for_genes <- match(temp$geneName,new_voom$geneName) #to get row ids of genes which are passed from voom translation and associated with repeats
  prepLMdata_genes<-na.omit(new_voom[v_ids_for_genes,])
  prepLMdata_genes<-unique(prepLMdata_genes)

  voom_data<-cbind(vall$geneName,v$E)
  writingResultOfVoom(voom_data,prefix)

  l<-list() # for last linear model list
  ncov<-0
  unique.repeats<-unique(prepLMdata_repeats$geneName)
  for (r in 1:length(unique.repeats)) {
    hit1<-prepLMdata_genes$geneName == as.character(unique.repeats[r])
    s1<-prepLMdata_genes[hit1,2:ncol(vall)]  # for gene counts
    hit2<-prepLMdata_repeats$geneName == as.character(unique.repeats[r])
    s2<-prepLMdata_repeats[hit2,] # for repeat counts
    if(nrow(s1)==1 & nrow(s2)>0){
      goalsMenu <- as.character(s2$repeatName)
      output <- as.data.frame(matrix(rep(0, 1 + length(goalsMenu)), nrow=1))
      names(output) <- c(unique(as.character(s2$geneName)), as.vector(s2$repeatName))
      output[1:length(s1),1]<-as.numeric(s1)
      for (var in 1:(ncol(output)-1)){
        output[1:length(s2[var,4:ncol(s2)]),var+1]<- as.numeric(s2[var,4:ncol(s2)])
      }
      if(!is.null(covariates)){
        ncov<-ncol(covariates)
        if(nrow(covariates)==(ncol(count.matrix)-1)){
          output<-cbind(output,covariates)
        }
      }
      colnames(output)[1]<-gsub("-","_",colnames(output)[1])
      if(r>1){
        formula <- paste(as.character(colnames(output)[1]),"~.",sep = "")
        lm.out<-lm(formula = formula , data=output)
        l<-list.append(l,lm.out)
      }else if (r==1){
        formula <- paste(as.character(colnames(output)[1]),"~.",sep = "")
        lm.out<-lm(formula = formula , data=output)
        l<-list.append(list(lm.out))
      }
    }

  }
  writingResultOfLM(l,ncov,prefix)
  return(l)
  }


writingResultOfLM<-function(lm_list,ncov,prefix){
  y<-data.frame()
  y<-as.data.frame(matrix(ncol=6,nrow=length(lm_list)))
  names(y) <- c("GeneName","RepeatName" , "r.squared" , "adjusted-r.squared" , "model-p.value", "individual-p.vals")
  id<-1
  for (list in lm_list) {
    y$GeneName[id]<-colnames(list$model)[1]
    y$RepeatName[id]<-paste(colnames(list$model)[2:(ncol(list$model)-ncov)],collapse = " ")
    y$r.squared[id]<-summary(list)$r.squared
    y$`adjusted-r.squared`[id]<-summary(list)$adj.r.squared
    y$`model-p.value`[id]<-lmp(list)
    n<-summary(list)$coefficients[,4]
    na<-names(list$coefficients[-1])
    if(is.logical(setdiff(na,names(n)))==0){
      y$`individual-p.vals`[id]<-paste(paste(names(n)[-1], n[-1], sep = " : ", collapse = " // "))
    }else{
      naList<-paste(setdiff(na,names(n)),": NA")
      y$`individual-p.vals`[id]<-paste(paste(names(n)[-1], n[-1], sep = " : ", collapse = " // "), paste(naList,collapse = " // "),sep = " // ")
    }
    id<-id+1
    if(lmp(list) < 0.05 & !is.na(lmp(list))){
      dir<-paste0(getwd(),gsub(" ","", paste("/",prefix,"-output/", colnames(list$model)[1])),collapse="")
      dir.create(dir,recursive = T)
      filedir<-paste0(dir,"/lmfitOTONE%1d.png",collapse = "")
      png(filedir, width=720, height=720, pointsize=16)
      plot(list)
      dev.off()
    }
  }
  write.table(y, file=paste(prefix,"-lm-results.tsv",collapse = ""), quote=F, sep="\t", row.names=F, col.names=T)
}

lmp <- function (modelobject){
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


writingResultOfVoom<-function(v,prefix){
  write.table(v, paste(prefix,"-cpm-values.tsv",collapse = ""), quote=F, sep="\t",row.names = F,col.names = T)
}
