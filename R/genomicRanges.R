
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
  gr <- with(df, GenomicRanges::GRanges(chr, IRanges(start, end), strand = strand))
  values(gr)<-df[,5:length(df)]
  values(gr)<-cbind(values(gr),df[2:3])
  return(gr)
}

makeGrObj_Unstrand <- function(df){  #within strandness
  gr <- with(df, GenomicRanges::GRanges(chr, IRanges(start, end), strand = "*"))
  values(gr)<-df[,5:length(df)]
  values(gr)<-cbind(values(gr),df[2:3])
  return(gr)
}

toFindOverlaps<-function(gr_repeats,gr_genome,strand){  #to get overlap in general function
  # Find overlaps
  names(values(gr_genome))[3:4]<-c("gStart","gEnd")
  names(values(gr_repeats))[4:5]<-c("rStart","rEnd")
  if(strand=="same"){
    m <- GenomicRanges::findOverlaps(gr_genome, gr_repeats)
    gr_genome.matched <- gr_genome[queryHits(m)]
    # Add the metadata from gr2
    mcols(gr_genome.matched) <- cbind.data.frame(
      mcols(gr_genome.matched),
      mcols(gr_repeats[subjectHits(m)]))

    return(gr_genome.matched)
  }else if (strand=="strandness"){
    m <- GenomicRanges::findOverlaps(gr_genome, gr_repeats,ignore.strand=TRUE)
    gr_genome.matched <- gr_genome[queryHits(m)]
    # Add the metadata from gr2
    mcols(gr_genome.matched) <- cbind.data.frame(
      mcols(gr_genome.matched),
      mcols(gr_repeats[subjectHits(m)]))

    return(gr_genome.matched)
  }


}

