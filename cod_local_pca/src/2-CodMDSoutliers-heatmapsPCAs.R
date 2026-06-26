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
  colnames(G) <- round(map$physical.pos/1000000,3)
  
  map <- my_snps$map
  pops <- my_snps$fam[-remove,]
  dim(pops)
  head(pops)
  
  rm(my_snps)

## Read in pop data and label rows of G
  popdata <- read.table("../../data/1-Pops.txt", header=TRUE)
  popdata$pop <- paste(popdata$region, popdata$ecotype, sep="-")
  popdata$colorEcotypes <- c(rep("lightblue", 5),
                             "darkviolet","cornflowerblue",
                             "darkviolet", "cornflowerblue")
  popdata
  pops <- merge(pops, popdata, by.x="family.ID", by.y="pop.names", all.x=TRUE)
  dim(pops)
  head(pops)
  pops$pop.sampleID <- paste(pops$pop, pops$sample.ID, sep="-")
  head(pops)
  rownames(G) <- pops$pop.sampleID # label G row names

# Reorganize G according to populations
  my.order <- order(pops$pop.sampleID)
  pops <- pops[my.order,] 
  
  G <- G[my.order,]
  
  identical(rownames(G), pops$pop.sampleID) # sanity check

# Read in the unique window table ####
  window_size=100
  windows <- read.table(paste0("../data_outputs/local_pca_results_UniquePCAwindows-",
                               window_size,"SNPwindow.txt"), header=TRUE)
  head(windows)

#### Read in the unique SV table ####
  SVs <- read.csv("../../outputs/INV-dataset_ConsolidationMethod1.5_unique-filtered.csv")
  head(SVs)
  SVs$invSizeMin[SVs$invSizeMin < 0] <- 0
  
  ## Add color by number of individuals
  color_endpoints <- c(grey(0.95), "#e2925d")
  SVs$color <- colByValue(SVs$n_ind_startWindows, colors = c(color_endpoints))
  
  summary(SVs$n_ind_startWindows)
  color_gradient <- t(colorRampPalette(color_endpoints)(100))

  pdf("../figures/ColorRampLegend.pdf", width=6, height=4)
  plot(x=1,y=1, col="white", xlim=c(0,100), ylim=c(0,8))
  rasterImage(color_gradient, 
              xleft = 10, ybottom = 4, 
              xright = 90, ytop = 5)
  rect(xleft = 10, ybottom = 4, xright = 90, ytop = 5, border = "black", lwd = 1)
  # 5. Add tick mark labels underneath the bar
  # 'adj' controls horizontal alignment (0 = left, 0.5 = center, 1 = right)
  vals_min=15
  vals_max=294
  text(x = 10, y = 3.5, labels = paste0(vals_min), adj = 0, cex = 0.9)
  text(x = 50, y = 3.5, labels = floor(vals_max-vals_min)/2, adj = 0.5, cex = 0.9)
  text(x = 90, y = 3.5, labels = paste0(vals_max), adj = 1, cex = 0.9)
  text(x = 50, y = 5.5, labels = "Number of individuals", font = 2, cex = 1.1)
  dev.off()
  

 #### plot inversion sizes ####
  pdf("../figures/InvSizeFreq.pdf", width=6, height=5)
  par(mar=c(4,4,1,1), mfrow=c(1,1), oma=c(0,0,0,0))
  hist(log10(SVs$invSizeMax), col="#e2925d", breaks=(seq(3, 8,by=0.1)), 
       xlab= "Log10(Inversion size in bases)", main="")
  hist(log10(windows$end-windows$start), breaks=(seq(3, 8,by=0.1)), add=TRUE, col="#5dade2")
  legend("topright",
         fill=c("#e2925d", "#5dade2"), legend=c("SV-based", "PCA-based"))
  
  par(mar=c(4,4,1,1), mfrow=c(1,1), oma=c(0,0,0,0))
  hist((SVs$invSizeMax)/10^6, col="cornflowerblue", breaks=(seq(0, 43,by=0.1)), 
       xlab= "Log10(Inversion size in bases)", main="")
  hist((windows$end-windows$start)/10^6, breaks=(seq(0, 43,by=0.1)), add=TRUE)
  legend("topright",
         fill=c("cornflowerblue", "grey"), legend=c("SV-based", "PCA-based"))
  dev.off()
  
  summary((windows$end-windows$start))
  summary(SVs$invSizeMax)


#### Plot windows ####
for (i in 1:nrow(windows)){
  print(i)
  subSNPs <- which(map$chromosome==windows$chromosome[i] & 
                     map$physical.pos > windows$start[i] &
                     map$physical.pos < windows$end[i]
  )
  nsnps <- length(subSNPs)
  Gsub <- G[,subSNPs]
  dim(Gsub)
  mainName <- paste0("WindowID: ", windows$PCAwindowID[i], "; ", 
                     windows$chromosome[i], " (LG ",windows$chr_num[i],"); ",
                     nsnps, " SNPs")
  pdf(paste0("../results-heatmapsPCAsByMDSOutlier-100windows/",
             windows$chromosome[i],"_Chrom",windows$chr_num[i],"_WindowID_",windows$PCAwindowID[i],".pdf"  ),
      width=6, height=6)
  
  # heatmaps can handle about 1000 snps
    if (ncol(Gsub)>3000){ Gplot <- Gsub[,sort(sample(1:ncol(Gsub), 3000))]}else{
      Gplot <- Gsub
    }
  
  ### Heatmap ####
  heatmap(Gplot, Colv = NA, scale = "none",
          main = mainName, RowSideColors = pops$colorEcotypes, cexRow=0.2
          )
  heatmap(Gplot, Colv = NA, Rowv = NA, scale = "none",
                   main = mainName, RowSideColors = pops$colorEcotypes,
          cexRow=0.2
  )
  
  #### PCA ####
  pcasub <- pca(Gsub, method = "ppca", nPcs = 2)
    #perform probabilitistic PCA for missing data
  scores <- scores(pcasub)
  head(scores)
  plot(scores, main=mainName, col=pops$colorEcotypes, 
       pch=as.numeric(factor(pops$colorEcotypes)))
  
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
  kpAddBaseNumbers(kp.plot, tick.dist=(windows$end[i]-windows$start[i])/5, 
                   tick.len=7.5, add.units=TRUE, units="Mb", cex=1, 
                   minor.tick.dist=(windows$end[i]-windows$start[i])/10)
  
  dev.off()
  
}



#### prep for karyoploter ####
#Create custom genome
gadmor3_df <- data.frame(chr=c("NC_044048.1", "NC_044049.1", "NC_044050.1", "NC_044051.1", "NC_044052.1", "NC_044053.1", "NC_044054.1", "NC_044055.1", "NC_044056.1", "NC_044057.1", "NC_044058.1", "NC_044059.1", "NC_044060.1", "NC_044061.1", "NC_044062.1", "NC_044063.1", "NC_044064.1", "NC_044065.1", "NC_044066.1", "NC_044067.1", "NC_044068.1", "NC_044069.1", "NC_044070.1"), start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), end=c(30875876, 28732775, 30954429, 43798135, 25300426, 27762770, 34137969, 29710654, 26487948, 27234273, 30713045, 30948897, 28829685, 29586942, 28657694, 34794352, 21723002, 24902675, 22015597, 24843429, 22358821, 23744039, 25242006))
head(gadmor3_df)
gadmor3genome <- toGRanges(gadmor3_df)
head(gadmor3genome)


### Karyotyper code for a specific region


my_plots <- function(chr, zoom_start, zoom_end){
  
  pdf(paste0("../results-heatmap-otherRanges/",
             windows$chromosome[i],"_Chrom",windows$chr_num[i],
             "start",zoom_start/10^6,"-end", zoom_end/10^6,"Mb.pdf"), width=6, height=6)
  whichSVs <- which(SVs$Chrom==chr #& 
                    #  (SVs$window_StartLowerCI < windows$end[i] |
                    #  SVs$window_EndUpperCI > windows$start[i])
  )
  zoom_region <- GRanges(seqnames = chr, 
                         ranges = IRanges(start = zoom_start, end = zoom_end))
  df <- SVs[whichSVs, c("Chrom", "window_StartLowerCI", "window_EndUpperCI", "color")]
  localSVs <- toGRanges(df)
  kp.chrom <- plotKaryotype(plot.type=1, genome=gadmor3genome, chromosome=chr, labels.plotter=NULL,
                            zoom=zoom_region)
  kpPlotRegions(kp.chrom, data=localSVs, col=localSVs$color, border=darker("#5dade2"), r0=0, r1=1)
  kpAddBaseNumbers(kp.chrom, tick.dist=(zoom_end-zoom_start)/5, tick.len=7.5, add.units=TRUE, 
                   units="Mb", minor.tick.dist=(zoom_end-zoom_start)/10, minor.tick.len=5, cex=1)
  

  
  
  #heatmap and pca
  subSNPs <- which(map$chromosome==chr  & 
                     map$physical.pos > zoom_start &
                     map$physical.pos < zoom_end
  )
  nsnps <- length(subSNPs)
  Gsub <- G[,subSNPs]
  # heatmaps can handle about 1000 snps
  if (ncol(Gsub)>3000){ Gplot <- Gsub[,sort(sample(1:ncol(Gsub), 3000))]}else{
    Gplot <- Gsub
  }
  
  ### Heatmap ####
  heatmap(Gplot, Colv = NA, scale = "none",
          main = mainName, RowSideColors = pops$colorEcotypes, cexRow=0.2
  )
  heatmap(Gplot, Colv = NA, Rowv = NA, scale = "none",
          main = mainName, RowSideColors = pops$colorEcotypes,
          cexRow=0.2
  )
  #### PCA ####
  pcasub <- pca(Gsub, method = "ppca", nPcs = 2)
  #perform probabilitistic PCA for missing data
  scores <- scores(pcasub)
  head(scores)
  plot(scores, col=pops$colorEcotypes, 
       pch=as.numeric(factor(pops$colorEcotypes)))
  dev.off()
 }


my_plots(chr = "NC_044051.1", zoom_start <- 17*10^6, zoom_end <- 19*10^6)
