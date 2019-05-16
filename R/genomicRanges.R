library(dplyr)
require(GenomicRanges)

#adfalsfmlsafa

getUpstream <- function(df,length, isWithGeneBody){ # to get genome interval of up stream of given genome

  dftemp <- df
  dfnew <- df

  indexP <- dftemp$strand == "+"
  indexN <- dftemp$strand == "-"

  dfnew$start[indexP] <- c(dftemp$start[indexP] - length)
  dfnew$end[indexN] <- c(dftemp$end[indexN] + length)

  if(isWithGeneBody == FALSE)
  {
    dfnew$end[indexP] <- c(dftemp$start[indexP])
    dfnew$start[indexN] <- c(dftemp$end[indexN] )
  }

  return(dfnew)
}

getDownstream <- function(df,length,isWithGeneBody){ # to get genome interval of down stream of given genome
  dftemp <- df
  dfnew <- df
  indexP <- dftemp$strand == "+"
  indexN <- dftemp$strand == "-"
  dfnew$start[indexN] <- dftemp$start[indexN] - length
  dfnew$end[indexP] <- dftemp$end[indexP] + length
  if(isWithGeneBody == FALSE)
  {
    dfnew$end[indexN] <- dftemp$start[indexN]
    dfnew$start[indexP] <- dftemp$end[indexP]
  }
  return(dfnew)
}

getDownAndUpStream <- function(df,len_up,len_down){ # to get genome interval of up and down stream of given genome with its
  dftemp <- df
  dfnew <- df

  indexP <- dftemp$strand == "+"
  indexN <- dftemp$strand == "-"

  dfnew$start[indexP] <- dftemp$start[indexP]-len_up
  dfnew$end[indexP] <- dftemp$end[indexP]+len_down

  dfnew$start[indexN] <- dftemp$start[indexN]-len_down
  dfnew$end[indexN] <- dftemp$end[indexN]+len_up

  return(dfnew)
}


makeGRangeObj <- function(df){  # to make grange object to analysis genome arithmetics well
  library(GenomicRanges)
    gr <- with(df, GRanges(chr, IRanges(start, end), strand = strand))
    values(gr)<-df[,5:length(df)]
  return(gr)
}

makeGrObj_Unstrand <- function(df){  #within strandness
    gr <- with(df, GRanges(chr, IRanges(start, end), strand = "*"))
    values(gr)<-df[,5:length(df)]
  return(gr)
}


toFindOverlaps<-function(gr_repeats,gr_genome){  #to get overlap in general function
  hits<-findOverlaps(gr_repeats, gr_genome)
  overlaps<-pintersect(gr_repeats[queryHits(hits)], gr_genome[subjectHits(hits)])
  return(overlaps)
}
