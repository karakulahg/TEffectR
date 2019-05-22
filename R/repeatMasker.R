

library(googledrive)
readRepeatMasker<-function(assembly){
  if(assembly=="hg19" || assembly=="Grch37"){
    dt<-drive_download(
      as_id("14-CbLlCgnYRr_GX6jiu-Ixn-5MziNcnQ"), type = ".csv", overwrite = TRUE) #https://drive.google.com/open?id=14-CbLlCgnYRr_GX6jiu-Ixn-5MziNcnQ
    # dt<-read.csv(dt$local_path,sep = "\t")
    return(dt)
  }else if(assembly=="hg38" || assembly=="Grch38"){
    dt<-drive_download(
      as_id("1PXVXY5_xC892QeuwXLoOscZnrzHeZG2B"), type = ".csv", overwrite = TRUE) #https://drive.google.com/open?id=1PXVXY5_xC892QeuwXLoOscZnrzHeZG2B
    # dt<-read.csv(dt$local_path,sep = "\t")
    return(dt)
  }else if(y=="mm9" || y =="Grcm37"){
    dt<-drive_download(
      as_id("1Lr977gf-sRnFQYxFA4TwXd9MBtnvQbjH"), type = ".csv", overwrite = TRUE) #https://drive.google.com/open?id=1Lr977gf-sRnFQYxFA4TwXd9MBtnvQbjH
    # dt<-read.csv(dt$local_path,sep = "\t")
    return(dt)
  }else if(y=="mm10" || y =="Grcm38"){
    dt<-drive_download(
      as_id("1qun2LN_5CmVwYFQ1PZww2GYFSuPog77F"), type = ".csv", overwrite = TRUE) #https://drive.google.com/open?id=1qun2LN_5CmVwYFQ1PZww2GYFSuPog77F
    # dt<-read.csv(dt$local_path,sep = "\t")
    return(dt)
  }else{
    print("not found assembly")
    return(NULL)
  }
}
