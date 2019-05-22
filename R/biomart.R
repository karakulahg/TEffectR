library(stringr)
library(biomaRt)
#for humans
#listDatasets(ensembl)


filterTranscriptID<-function(transcript_ids,assembly){
  if(length(transcript_ids)!=0){
    ensembl = returnEnsembl(assembly)
    if(!is.null(ensembl)){
      data <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                   "chromosome_name", "start_position", "end_position",
                                   "strand"),

                    filters = c("ensembl_transcript_id"), #listFilters(ensembl)
                    values = transcript_ids,
                    mart=ensembl,
                    verbose = TRUE)

      names(data)<-c("geneID","geneName","chr","start","end","strand")
      data <- data[c(3,4,5,6,1,2)]

      if(length(data)>0){
        hit<-stringr::str_detect(data$strand,"-1")
        data$strand[hit]<-"-"

        hit<-stringr::str_detect(data$strand,"1")
        data$strand[hit]<-"+"

      }

      return(data)
    }

  }else return(NULL)


}

filterTranscriptID_V<-function(transcript_ids,assembly){
  if(length(transcript_ids)!=0){
    ensembl = returnEnsembl(assembly)
    if(!is.null(ensembl)){
      data <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                   "chromosome_name", "start_position", "end_position",
                                   "strand"),

                    filters = c("ensembl_transcript_id_version"), #listFilters(ensembl)
                    values = transcript_ids,
                    mart=ensembl,
                    verbose = TRUE)

      names(data)<-c("geneID","geneName","chr","start","end","strand")
      data <- data[c(3,4,5,6,1,2)]

      if(length(data)>0){
        hit<-stringr::str_detect(data$strand,"-1")
        data$strand[hit]<-"-"

        hit<-stringr::str_detect(data$strand,"1")
        data$strand[hit]<-"+"
      }

      return(data)
    }

  }else return(NULL)


}


filterGeneName<-function(genes,assembly){
  if(length(genes)!=0){
    ensembl = returnEnsembl(assembly)
    gene_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = c("hgnc_symbol"), values = genes, mart = ensembl, verbose = T)
    if(!is.null(ensembl) & !is.null(gene_ids)){
      data <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                   "chromosome_name", "start_position", "end_position",
                                   "strand"),

                    filters = c("ensembl_gene_id"), #listFilters(ensembl)
                    values = gene_ids$ensembl_gene_id,
                    mart=ensembl,
                    verbose = TRUE)

      names(data)<-c("geneID","geneName","chr","start","end","strand")
      data <- data[c(3,4,5,6,1,2)]

      if(length(data)>0){
          hit<-stringr::str_detect(data$strand,"-1")
          data$strand[hit]<-"-"

          hit<-stringr::str_detect(data$strand,"1")
          data$strand[hit]<-"+"
      }


      return(data)
    }

  }else return(NULL)

}


filterGeneID_V<-function(gene_ids,assembly){
  if(length(gene_ids)!=0){
    ensembl = returnEnsembl(assembly)
    if(!is.null(ensembl) & !is.null(gene_ids)){
      data <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                   "chromosome_name", "start_position", "end_position",
                                   "strand"),

                    filters = c("ensembl_gene_id_version"), # listFilters(ensembl)
                    values = gene_ids,
                    mart=ensembl,
                    verbose = TRUE)

      names(data)<-c("geneID","geneName","chr","start","end","strand")
      data <- data[c(3,4,5,6,1,2)]

      if(length(data)>0){
        hit<-stringr::str_detect(data$strand,"-1")
        data$strand[hit]<-"-"

        hit<-stringr::str_detect(data$strand,"1")
        data$strand[hit]<-"+"
      }

      return(data)
    }

  }else return(NULL)

}


filterGeneID<-function(gene_ids,assembly){
  if(length(gene_ids)!=0){
    ensembl = returnEnsembl(assembly)
    if(!is.null(ensembl) & !is.null(gene_ids)){
      data <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                   "chromosome_name", "start_position", "end_position",
                                   "strand"),

                    filters = c("ensembl_gene_id"), # listFilters(ensembl)
                    values = gene_ids,
                    mart=ensembl,
                    verbose = TRUE)

      names(data)<-c("geneID","geneName","chr","start","end","strand")
      data <- data[c(3,4,5,6,1,2)]

      if(length(data)>0){
        hit<-stringr::str_detect(data$strand,"-1")
        data$strand[hit]<-"-"

        hit<-stringr::str_detect(data$strand,"1")
        data$strand[hit]<-"+"
      }

      return(data)
    }

  }else return(NULL)

}




returnEnsembl<-function(assembly){library(biomaRt)
  if(assembly=="hg19" || assembly=="Grch37"){
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = 37, verbose = TRUE)
    return(ensembl)
  }else if(assembly=="hg38" || assembly=="Grch38"){
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = 38, verbose = TRUE)
    return(ensembl)
  }else if(y=="mm9" || y =="Grcm37"){
    ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", verbose = TRUE)
    return(ensembl)
  }else if(y=="mm10" || y =="Grcm38"){
    ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", verbose = TRUE)
    return(ensembl)
  }else{
    print("not found assembly")
    return(NULL)}
}





