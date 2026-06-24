library(rlang)
#devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
library(bigstatsr)
library(bigsnpr)
library(Matrix)
library(dplyr)
library(LEA)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")
library(pcaMethods)

setwd("~/Desktop/cod_local_pca/src")


my_snps <- snp_attach("../data/merged.f.99ind.MAF05.rds")
head(my_snps)
G <- my_snps$genotypes
G <- G[-17,]
dim(G)
str(G)
head(G[,1:5])
map <- my_snps$map
pops <- my_snps$fam[-17,]
rownames(G) <- pops$family.ID
colnames(G) <- round(map$physical.pos/1000000,3)
rm(my_snps)

windows <- read.table("../data_outputs/local_pca_results_AllPCAwindows.txt", header=TRUE)
head(windows)

for (i in 1:nrow(windows)){
  print(i)
  subSNPs <- which(map$chromosome==windows$chromosome[i] & 
                     map$physical.pos > windows$start[i] &
                     map$physical.pos < windows$end[i]
  )
  length(subSNPs)
  Gsub <- G[,subSNPs]
  dim(Gsub)
  mainName <- paste0("Window ID ", windows$PCAwindowID[i], " ", 
                     windows$chromosome[i], " (",windows$chr_num[i],")")
  pdf(paste0("../results-heatmapsPCAsByMDSOutlier/",
             windows$chromosome[i],"_Chrom",windows$chr_num[i],"_WindowID_",windows$PCAwindowID[i],".pdf"  ),
      width=6, height=6)
  # heatmaps can handle about 1000 snps
    if (ncol(Gsub)>3000){ Gplot <- Gsub[,sort(sample(1:ncol(Gsub), 3000))]}else{
      Gplot <- Gsub
    }
  heatmap(Gplot, Colv = NA, scale = "none",
          main = mainName
          )
  pcasub <- pca(Gsub, method = "ppca", nPcs = 2)
    #perform probabilitistic PCA for missing data
  scores <- scores(pcasub)
  head(scores)
  plot(scores, main=mainName)
  dev.off()
}

### Plot histogram of PCA window sizes

subSNPs <- which(map$chromosome==windows$chromosome[1] & 
                   map$physical.pos > 5.4*10^6 &
                   map$physical.pos < 5.6*10^6
)
length(subSNPs)
Gsub <- G[,subSNPs]
dim(Gsub)
heatmap(Gsub, Colv = NA, scale = "none",
        main = mainName
)

hist(log10(windows$end-windows$start), breaks=(seq(0, 8,by=0.1)))
