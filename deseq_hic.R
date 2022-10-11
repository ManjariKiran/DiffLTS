#Loading all the ;ibraries
library(strawr)
library(tidyverse)
library(dbscan)
library(InteractionSet)
library(raster)
library(devtools)
#install_github('dozmorovlab/HiCcompare', build_vignettes = TRUE)
library(HiCcompare)
library(hictoolsr)
library(terra)

setwd("/Users/phanstiel2/MK/Work")


#fs_loops <- read.delim("/Users/phanstiel2/MK/Work/FS_5kbLoops.txt")

## Merge loops and convert to GInteractions

loops <- mergeBedpe(bedpeFiles = c("WT_5kbLoops.txt", "FS_5kbLoops.txt"), res = 10e3) |> as_ginteractions()

head(loops)

setwd("/Users/phanstiel2/MK/Work")
##Extracting Hi-C counts
hicFiles <- c("GSM4259900_HEK_HiC_NUP_IDR_FS_A9_1_1_inter_30.hic","GSM4259901_HEK_HiC_NUP_IDR_FS_A9_1_2_inter_30.hic","GSM4259902_HEK_HiC_NUP_IDR_FS_A9_2_1_inter_30.hic","GSM4259903_HEK_HiC_NUP_IDR_FS_A9_2_2_inter_30.hic","GSM4259896_HEK_HiC_NUP_IDR_WT_A9_1_1_inter_30.hic","GSM4259897_HEK_HiC_NUP_IDR_WT_A9_1_2_inter_30.hic","GSM4259898_HEK_HiC_NUP_IDR_WT_A9_2_1_inter_30.hic","GSM4259899_HEK_HiC_NUP_IDR_WT_A9_2_2_inter_30.hic")
hicFiles

loopCounts <- extractCounts(bedpe = loops,
                            hic = hicFiles,
                            chroms = c(1:22, "X"),
                            res = 10e3,
                            norm = 'NONE',
                            matrix = 'observed')

library(InteractionSet)
## Simplify column names

colnames(mcols(loopCounts)) <- 
  gsub(pattern = "GSM.*_IDR_(WT|FS)_A9_(1|2)_(1|2)_.*", 
       replacement = "\\1_\\2_\\3",
       x = colnames(mcols(loopCounts)))
head(loopCounts)
library(DESeq2)
## Isolate count matrix
cnts <- 
  mcols(loopCounts)[grep("WT|FS", colnames(mcols(loopCounts)))] |>
  as.matrix()
head(cnts)

colData <- 
  do.call(rbind, strsplit(x = colnames(cnts), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "biorep", "techrep"))
colData

dds <- 
  DESeqDataSetFromMatrix(countData = round(cnts),
                         colData = colData,
                         design = ~techrep + biorep + condition)

dds <- DESeq(dds)
sizeFactors(dds)


## Run DEseq analysis
res <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_WT_vs_FS", type="apeglm")
summary(res) 

## Attach DESeq2 results
mcols(loopCounts) <- cbind(mcols(loopCounts), res)
## Separate WT/FS-specific loops
wtLoops <- loopCounts[loopCounts$padj <= 0.01 &
                        loopCounts$log2FoldChange > 0]
fsLoops <- loopCounts[loopCounts$padj <= 0.01 &
                        loopCounts$log2FoldChange < 0]

write.table(loopCounts,"WT_vs_FS_no_norm",quote=FALSE,sep="\t")
###For the normalized expected data

# make object

library(DESeq2)

samples <- c("GSM4259900_HEK_HiC_NUP_IDR_FS_A9_1_1_inter_30.hic","GSM4259901_HEK_HiC_NUP_IDR_FS_A9_1_2_inter_30.hic","GSM4259902_HEK_HiC_NUP_IDR_FS_A9_2_1_inter_30.hic","GSM4259903_HEK_HiC_NUP_IDR_FS_A9_2_2_inter_30.hic","GSM4259896_HEK_HiC_NUP_IDR_WT_A9_1_1_inter_30.hic","GSM4259897_HEK_HiC_NUP_IDR_WT_A9_1_2_inter_30.hic","GSM4259898_HEK_HiC_NUP_IDR_WT_A9_2_1_inter_30.hic","GSM4259899_HEK_HiC_NUP_IDR_WT_A9_2_2_inter_30.hic")

mh_index <- function(buffer, loop, inner){
  m=(buffer*2)+1
  center <- buffer+1
  M <- matrix(data=NA,nrow=m,ncol=m)
  for(j in 1:m){
    l=buffer+j
    for(i in 1:m){
      k=(m+1)-j
      if((i <= (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-k-i
      }
      if((i <= (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-j-i
      }
      if((i > (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-i-j
      }
      if((i > (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-i-k
      }
    }
  }
  inner_val <- which(M == inner)
  new <- loop[inner_val]
  return(new)
}


counts <- cnts
expected=data.frame()
normalized=data.frame()

for(i in 1:length(loops)){
  for(j in 1:length(hicFiles)){
 # l <- calcApa(bedpe = loops[j], hic = hicFiles[i], norm = "NONE", res = 10e3, buffer = 10, 
 #              filter = FALSE)
  loe <- calcApa(bedpe = loops[i], hic = hicFiles[j], norm = "NONE", res = 10e3, buffer = 10, 
                filter = FALSE, matrix = 'oe')
  for(k in 1:3){
    inner = k + 1
    normalized[i,k] <- median(mh_index(buffer = 10, loop = loe, inner))
    #normalized[j,i+1]<- (median(mh_index(buffer = 10, loop = l, inner = i))+1)/(median(mh_index(buffer = 10, loop = le, inner = i))+1)
  }
  expected[i,j] <- round(rowMeans(normalized[i,]),3)
  }
}
dim(expected)

colnames(expected) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2","WT_1_1","WT_1_2","WT_2_1","WT_2_2")
head(cnts)
head(expected)
loopobject = list()
loopobject[[1]] = counts
loopobject[[2]] = expected
#loopobject[[5]] = distances
#loopobject[[6]] = norm_depth
#loopobject[[7]] = norm_dist
loopobject[[3]] = samples


norm_MH = loopobject[[2]]
norm_MH
# make DDS object
dds <- 
  DESeqDataSetFromMatrix(countData = round(cnts),
                         colData = colData,
                         design = ~techrep + biorep + condition)


head(norm_MH)
# Add the normalization factors
# center values around 1

norm_MH = as.matrix(norm_MH) / exp(rowMeans(log(as.matrix(norm_MH))))
norm_MH[!is.finite(norm_MH)] <- 1
#norm_MH[!complete.cases(norm_MH),]
dim((norm_MH))
# add normalization values
normalizationFactors(dds) <- as.matrix(norm_MH)

# perform enrichments
dds <- DESeq(dds,betaPrior=FALSE, fitType="local")
res <- results(dds)

res

dds <- DESeq(dds)
sizeFactors(dds)


## Run DEseq analysis
res <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_WT_vs_FS", type="apeglm")
summary(res) 

## Attach DESeq2 results
mcols(loopCountse) <- cbind(mcols(loopCountse), res)
## Separate WT/FS-specific loops
wtLoops <- loopCountse[loopCountse$padj <= 0.01 &
                         loopCounts$log2FoldChange > 0]
fsLoops <- loopCountse[loopCountse$padj <= 0.01 &
                         loopCountse$log2FoldChange < 0]


write.table(loopCounts,"WT_vs_FS_1-3",quote=FALSE,sep="\t")





