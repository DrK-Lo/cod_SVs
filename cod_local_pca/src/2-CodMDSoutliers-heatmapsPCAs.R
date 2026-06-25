# Uncomment to install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pcaMethods")
# BiocManager::install("karyoploteR")
# BiocManager::install("BSgenome")

library(rlang)
#devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
library(bigstatsr)
library(bigsnpr)
library(Matrix)
library(dplyr)
library(LEA)
library(pcaMethods)
library(karyoploteR)
library(BSgenome)
library(ggplot2)

setwd("~/Documents/GitHub/cod_SVs/cod_local_pca/src")


my_snps <- snp_attach("~/Desktop/codGenotypeData/merged.f.99ind.MAF05.rds")
  # this file is too large for the GitHub

remove = c(17,19) # see last scripts for removal of these individuals
head(my_snps)
G <- my_snps$genotypes
G <- G[-remove,]
dim(G)
str(G)
head(G[,1:5])
map <- my_snps$map
pops <- my_snps$fam[-remove,]
rownames(G) <- pops$family.ID
colnames(G) <- round(map$physical.pos/1000000,3)

map <- my_snps$map
pops <- my_snps$fam[-remove,]

rm(my_snps)

# Read in the unique window table
window_size=100
windows <- read.table(paste0("../data_outputs/local_pca_results_UniquePCAwindows-",
                             window_size,"SNPwindow.txt"), header=TRUE)
head(windows)

#### Read in the unique SV table ####
SVs <- read.csv("../../outputs/INV-dataset_ConsolidationMethod1.5_unique-filtered.csv")
head(SVs)
SVs$invSizeMin[SVs$invSizeMin < 0] <- 0


#### plot inversion sizes ####
pdf("../figures/InvSizeFreq.pdf", width=6, height=5)
par(mar=c(4,4,1,1), mfrow=c(1,1), oma=c(0,0,0,0))
hist(log10(SVs$invSizeMax), col="cornflowerblue", breaks=(seq(3, 8,by=0.1)), 
     xlab= "Log10(Inversion size in bases)", main="")
hist(log10(windows$end-windows$start), breaks=(seq(3, 8,by=0.1)), add=TRUE)
legend("topright",
       fill=c("cornflowerblue", "grey"), legend=c("SV-based", "PCA-based"))

par(mar=c(4,4,1,1), mfrow=c(1,1), oma=c(0,0,0,0))
hist((SVs$invSizeMax)/10^6, col="cornflowerblue", breaks=(seq(0, 43,by=0.1)), 
     xlab= "Log10(Inversion size in bases)", main="")
hist((windows$end-windows$start)/10^6, breaks=(seq(0, 43,by=0.1)), add=TRUE)
legend("topright",
       fill=c("cornflowerblue", "grey"), legend=c("SV-based", "PCA-based"))
dev.off()

summary((windows$end-windows$start))
summary(SVs$invSizeMax)


#### prep for karyoploter ####
#Create custom genome
gadmor3_df <- data.frame(chr=c("NC_044048.1", "NC_044049.1", "NC_044050.1", "NC_044051.1", "NC_044052.1", "NC_044053.1", "NC_044054.1", "NC_044055.1", "NC_044056.1", "NC_044057.1", "NC_044058.1", "NC_044059.1", "NC_044060.1", "NC_044061.1", "NC_044062.1", "NC_044063.1", "NC_044064.1", "NC_044065.1", "NC_044066.1", "NC_044067.1", "NC_044068.1", "NC_044069.1", "NC_044070.1"), start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), end=c(30875876, 28732775, 30954429, 43798135, 25300426, 27762770, 34137969, 29710654, 26487948, 27234273, 30713045, 30948897, 28829685, 29586942, 28657694, 34794352, 21723002, 24902675, 22015597, 24843429, 22358821, 23744039, 25242006))
head(gadmor3_df)
gadmor3genome <- toGRanges(gadmor3_df)
head(gadmor3genome)

#### Plot windows ####
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
  pdf(paste0("../results-heatmapsPCAsByMDSOutlier-100windows/",
             windows$chromosome[i],"_Chrom",windows$chr_num[i],"_WindowID_",windows$PCAwindowID[i],".pdf"  ),
      width=6, height=6)
  # heatmaps can handle about 1000 snps
    if (ncol(Gsub)>3000){ Gplot <- Gsub[,sort(sample(1:ncol(Gsub), 3000))]}else{
      Gplot <- Gsub
    }
  
  ### Heatmap ####
  heatmap(Gplot, Colv = NA, scale = "none",
          main = mainName
          )
  
  #### PCA ####
  pcasub <- pca(Gsub, method = "ppca", nPcs = 2)
    #perform probabilitistic PCA for missing data
  scores <- scores(pcasub)
  head(scores)
  plot(scores, main=mainName)
  
  ### SVs ####
  whichSVs <- which(SVs$Chrom==windows$chromosome[i] #& 
                    #  (SVs$window_StartLowerCI < windows$end[i] |
                    #  SVs$window_EndUpperCI > windows$start[i])
                      )
  whichSVs
  df <- SVs[whichSVs, c("Chrom", "window_StartLowerCI", "window_EndUpperCI")]
  df$Chrom="."
  localSVs <- toGRanges(df)
  localGenome = toGRanges(data.frame(chr=".", start= windows$start[i], end = windows$end[i]))
  kp.plot <- plotKaryotype(genome=localGenome, chromosome=".")
  kpPlotRegions(kp.plot, data=localSVs, col="#5dade2", border=darker("#5dade2"), r0=0, r1=1)
  kpAddBaseNumbers(kp.plot, tick.dist=(windows$end[i]-windows$start[i])/4, 
                   tick.len=7.5, add.units=TRUE, units="Mb", cex=1, )
  
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


