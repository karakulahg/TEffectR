# # ###################################### <<<< 1 >>>> ###################################
# # ## To Read gene or transcript expression file
# x<-read.csv("~/Documents/G.Karakulah/BC/gene_count_matrix.csv", row.names = 1, header=T, stringsAsFactors = F)
# x<-x[1:250,]
# #
# # ######################################################################################
# #
# #
# # ###################################### <<<< 2 >>>> ###################################
# #

# library(stringr)
# library(biomaRt)
# library(biomartr)
# library(dplyr)
# library(Rsamtools)
# library(edgeR)
# library(rlist)
# library(limma)
#
#
# #
# # ## for gene annotation
# #
#    genes<-get_intervals(x = rownames(x), organism="hg38", ID.type = "ensembl_gene_id", URL="dec2014.archive.ensembl.org" ) #use organism parameter if you provide gene symbols
# #
# # ########################################################################################
# #
# # ###################################### <<<< 3 >>>> #####################################
# # #
# # # #for repeat annotation
# # # download repeatmasker which you will use it from website is "http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html"
# #
#    rm<-rm_format(filepath="~/Downloads/hg38.fa.out.gz")
# #
# # #########################################################################################
# #
# # ###################################### <<<< 4 >>>> ######################################
# # #
# # #
#    w<-get_overlaps(g = genes, r = rm, strand = "same", distance = 10000, repeat_class = "LINE")
# # #
# # # write.table(w, file="overlapped.bed", quote=F, sep="\t", row.names=F, col.names=F)
# # #
# # #
#
#    bamlist <- c("~/Documents/G.Karakulah/BC/uniqueBam/SRR5962198/SRR5962198_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962199/SRR5962199_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962200/SRR5962200_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962201/SRR5962201_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962202/SRR5962202_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962203/SRR5962203_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962204/SRR5962204_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962205/SRR5962205_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962206/SRR5962206_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962207/SRR5962207_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962208/SRR5962208_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962209/SRR5962209_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962210/SRR5962210_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962211/SRR5962211_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962212/SRR5962212_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962213/SRR5962213_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962214/SRR5962214_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962215/SRR5962215_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962216/SRR5962216_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962217/SRR5962217_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962218/SRR5962218_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962219/SRR5962219_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962220/SRR5962220_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962221/SRR5962221_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962222/SRR5962222_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962223/SRR5962223_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962224/SRR5962224_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962225/SRR5962225_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962226/SRR5962226_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962227/SRR5962227_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962228/SRR5962228_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962229/SRR5962229_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962230/SRR5962230_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962231/SRR5962231_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962232/SRR5962232_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962233/SRR5962233_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962234/SRR5962234_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962235/SRR5962235_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962236/SRR5962236_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962237/SRR5962237_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962238/SRR5962238_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962239/SRR5962239_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962240/SRR5962240_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962241/SRR5962241_unique_sorted.bam")
#
#
#    namelist<- c("SRR5962198","SRR5962199","SRR5962200","SRR5962201","SRR5962202","SRR5962203","SRR5962204","SRR5962205","SRR5962206","SRR5962207",
#                 "SRR5962208","SRR5962209","SRR5962210","SRR5962211","SRR5962212","SRR5962213","SRR5962214","SRR5962215","SRR5962216","SRR5962217",
#                 "SRR5962218","SRR5962219","SRR5962220","SRR5962221","SRR5962222","SRR5962223","SRR5962224","SRR5962225","SRR5962226","SRR5962227",
#                 "SRR5962228","SRR5962229","SRR5962230","SRR5962231","SRR5962232","SRR5962233","SRR5962234","SRR5962235","SRR5962236","SRR5962237",
#                 "SRR5962238","SRR5962239","SRR5962240","SRR5962241")
#
#
#    bamlist <- c("~/Documents/G.Karakulah/BC/uniqueBam/SRR5962198/SRR5962198_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962199/SRR5962199_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962200/SRR5962200_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962201/SRR5962201_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962202/SRR5962202_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962203/SRR5962203_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962204/SRR5962204_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962205/SRR5962205_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962206/SRR5962206_unique_sorted.bam",
#                 "~/Documents/G.Karakulah/BC/uniqueBam/SRR5962207/SRR5962207_unique_sorted.bam")
#
#    namelist<- c("SRR5962198","SRR5962199","SRR5962200","SRR5962201","SRR5962202","SRR5962203","SRR5962204","SRR5962205","SRR5962206","SRR5962207")
#
#    write.table(w, file="overlapped.csv", quote=F, sep="\t", row.names=F, col.names=T)
#
#    e<-rm_count(bamlist = bamlist, namelist, ranges=w) #for counting
#
#    summ<-co_summarise(e, namelist) #repeat.counts
#
#
#    write.table(e, file="counts-all-samples.csv", quote=F, sep="\t", row.names=F, col.names=T)
#
#    write.table(summ, file="summarized-all-samples.csv", quote=F, sep="\t", row.names=F, col.names=T)
#
#
#    genes.expr<-x
#    repeats.expr<-summ
#
#    merge_counts<-make_matrix(gene.annotation = genes, genes.expr = x[,1:10], repeats.expr = summ )
#
#    dge<-get_dge_transformation(count.matrix = count.matrix)
#
#    # for gene.counts :
#    # # df1<-genes[,5:6] ( # x<-read.csv("~/Downloads/gene_count_matrix.csv", row.names = 1, header=T, stringsAsFactors = F))
#    #
#    # # y<- data.frame(geneID = row.names(x), x)
#    # # df<-merge(df1,y,by="geneID")
#    # # View(df)
#    # # df<-df[,2:12] #gene.counts
#
#    s3<-c(rep("T",5), rep("N",5)) # to define a group
#    x <- data.frame("Data" =c(rep("T",5), rep("N",5)), "Age" = c(rep(13,2), rep(24,3), rep(32,4),23), "Sex" = c(rep("F",5),rep("M",5)))
#    glm_list<-apply_lm(dge.transformation = dge, repeat.counts = summ , group = s3)
#
#    write.table(last, file="matrix-all-samples.csv", quote=F, sep="\t", row.names=F, col.names=T)
#

# glm_list<-apply_lm(count.matrix,summ , x)

# glm_list[[5]]$model
# summary(glm_list[[4]])$r.squared #!
# summary(glm_list[[4]])$adj.r.squared #!
# summary(glm_list[[4]])$fstatistic
# summary(glm_list[[5]])$coefficients
# summary(glm_list[[5]])$coefficients[,4]
# summary(glm_list[[4]])
#
# lmp <- function (modelobject) {
#   if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
#   f <- summary(modelobject)$fstatistic
#   p <- pf(f[1],f[2],f[3],lower.tail=F)
#   attributes(p) <- NULL
#   return(p)
# }
#
# lmp(glm_list[[4]]) #!
#
#
#
# y<-data.frame()
# y<-as.data.frame(matrix(ncol=6,nrow=34))
# names(y) <- c("GeneName","Repeats-family/families" , "r.squared" , "adjusted-r.squared" , "model-p.value", "individual-p.vals")
# id<-1
# for (list in glm_list) {
#
#   y$GeneName[id]<-colnames(list$model)[1]
#   y$`Repeats-family/families`[id]<-paste(colnames(list$model)[2:(ncol(list$model)-ncol(covariates))],collapse = " ")
#   y$r.squared[id]<-summary(list)$r.squared
#   y$`adjusted-r.squared`[id]<-summary(list)$adj.r.squared
#   y$`model-p.value`[id]<-lmp(list)
#   if((ncol(list$model)-ncol(covariates))!=ncol(covariates)-1){
#     y$`indivudual-p.vals`[id]<-paste(paste(summary(glm_list[[id]])$coefficients[,4])[2:(ncol(list$model)-ncol(covariates))], collapse = " ")
#   }else{
#     y$`indivudual-p.vals`[id]<-"NA"
#   }
#
#   id<-id+1
#
# }
#
# write.table(y, file="results-lm.csv", quote=F, sep="\t", row.names=F, col.names=T)


