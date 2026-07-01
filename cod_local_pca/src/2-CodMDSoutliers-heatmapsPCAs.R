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
  print(my_snps$fam[remove,])
  pops <- my_snps$fam[-remove,]
  colnames(G) <- round(map$physical.pos/1000000,3)
  
  map <- my_snps$map
  pops <- my_snps$fam[-remove,]
  dim(pops)
  head(pops)
  
  rm(my_snps)

## Read in pop data and label rows of G ####
  popdata <- read.table("../../data/1-Pops.txt", header=TRUE)
  popdata$pop <- paste(popdata$region, popdata$ecotype, sep="-")
  popdata$colorEcotypes <- c(rep("lightblue", 5),
                             "darkviolet","cornflowerblue",
                             "darkviolet", "cornflowerblue")
  popdata
  write.table(popdata, file="../../data/1-PopsColors.txt")
  pops <- merge(pops, popdata, by.x="family.ID", by.y="pop.names", all.x=TRUE)
  dim(pops)
  head(pops)
  pops$pop.sampleID <- paste(pops$pop, pops$sample.ID, sep="-")
  head(pops)
  rownames(G) <- pops$pop.sampleID # label G row names

# Reorganize G according to populations ####
  my.order <- order(pops$pop.sampleID)
  pops <- pops[my.order,] 
  
  G <- G[my.order,]
  
  identical(rownames(G), pops$pop.sampleID) # sanity check

# Read in the unique window table ####
  window_size=100
  windows <- read.table(paste0("../data_outputs/local_pca_results_UniquePCAwindows-",
                               window_size,"SNPwindow.txt"), header=TRUE)
  head(windows)
  table(windows$chr_num)

#### Read in the unique SV table ####
    SVs <- read.csv("../../consolidatev2/data_outputs/INV-dataset-ConsensusID_unique.csv")
    head(SVs)
    #SVs$invSizeMin[SVs$invSizeMin < 0] <- 0 OLD OUTPUTS
    summary(SVs$n_PASS_ind_group)
    SVs <- SVs[which(SVs$pass2_ind>0),]
    dim(SVs)
    SVs <- SVs[order(SVs$chrom, SVs$start_minus_CI_group, SVs$end_plus_CI_group),]
  
  ## Add color by number of individuals and print legend ####
  color_endpoints <- c(grey(0.95), "#5dade2", "darkblue")
  #SVs$color <- colByValue(SVs$n_ind_invID, colors = c(color_endpoints)) # old code
  SVs$color <- colByValue(SVs$n_PASS_ind, 
                          colors = c(color_endpoints), min=0, max=294)
  
  
  color_gradient <- t(colorRampPalette(color_endpoints)(100))

  pdf("../figures/ColorRampLegend.pdf", width=6, height=4)
  plot(x=1,y=1, col="white", xlim=c(0,100), ylim=c(0,8))
  rasterImage(color_gradient, 
              xleft = 10, ybottom = 4, 
              xright = 90, ytop = 5)
  rect(xleft = 10, ybottom = 4, xright = 90, ytop = 5, border = "black", lwd = 1)
  # 5. Add tick mark labels underneath the bar
  # 'adj' controls horizontal alignment (0 = left, 0.5 = center, 1 = right)
  vals_min=0
  vals_max=294
  text(x = 10, y = 3.5, labels = paste0(vals_min), adj = 0, cex = 0.9)
  text(x = 50, y = 3.5, labels = floor(vals_max-vals_min)/2, adj = 0.5, cex = 0.9)
  text(x = 90, y = 3.5, labels = paste0(vals_max), adj = 1, cex = 0.9)
  text(x = 10, y = 5.5, labels = 0, adj = 0, cex = 0.9)
  text(x = 50, y = 5.5, labels = 0.5, adj = 0.5, cex = 0.9)
  text(x = 90, y = 5.5, labels = 1, adj = 1, cex = 0.9)
  text(x = 50, y = 2.8, labels = "Number of individuals", font = 2, cex = 1.1)
  text(x = 50, y = 6.2, labels = "Proportion of individuals", font = 2, cex = 1.1)
  dev.off()
  

### Read in individual SVs ####
  SVs_ind <- read.csv("../../consolidatev2/data_outputs/INV-dataset-ConsensusID-Ind-Filtered.csv")
  head(SVs_ind) 
  SVs_ind <- SVs_ind[which(SVs_ind$pass_ind),]
  head(SVs_ind)
  
  
 #### plot inversion sizes ####
  pdf("../figures/InvSizeFreq.pdf", width=6, height=5)
  par(mar=c(4,4,1,1), mfrow=c(1,1), oma=c(0,0,0,0))
  
  hist(log10(SVs$maxSize), col="#5dade2", breaks=(seq(2, 8,by=0.1)),
       xlab= "Log10(Inversion size in bases)", main="", ylim=c(0,55), cex=1.5)
  
  hist(log10(windows$end-windows$start), breaks=(seq(2, 8,by=0.1)), add=TRUE, 
       col=adjustcolor("#e2925d", alpha=0.8))
  legend("topright",
         fill=c("#5dade2", "#e2925d"), legend=c("Read-based", "PCA-based"),
         bty="n", cex=1.5)
  
  # SVs "#5dade2"
  # MDS "#e2925d"
  
  dev.off()
  
  summary((windows$end-windows$start))
  summary(SVs$maxSize)

#### prep for karyoploter ####
  #Create custom genome
  gadmor3_df <- data.frame(chr=c("NC_044048.1", "NC_044049.1", "NC_044050.1", "NC_044051.1", "NC_044052.1", "NC_044053.1", "NC_044054.1", "NC_044055.1", "NC_044056.1", "NC_044057.1", "NC_044058.1", "NC_044059.1", "NC_044060.1", "NC_044061.1", "NC_044062.1", "NC_044063.1", "NC_044064.1", "NC_044065.1", "NC_044066.1", "NC_044067.1", "NC_044068.1", "NC_044069.1", "NC_044070.1"), start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), end=c(30875876, 28732775, 30954429, 43798135, 25300426, 27762770, 34137969, 29710654, 26487948, 27234273, 30713045, 30948897, 28829685, 29586942, 28657694, 34794352, 21723002, 24902675, 22015597, 24843429, 22358821, 23744039, 25242006))
  head(gadmor3_df)
  gadmor3genome <- toGRanges(gadmor3_df)
  head(gadmor3genome)



### SV plotting function with zoom  ####
my_plots <- function(chr_num, zoom_start, zoom_end, folder, windowID=""){
  
  chr <- chromosomes[chr_num]
  mainName = chr
  pdf(paste0(folder,
             chr,"_Chrom-", chr_num, "_windowID-", windowID,
             "_start-", round(zoom_start/10^6,3),"_end-", round(zoom_end/10^6,3),"Mb.pdf"), width=6, height=6)

  
  chr = chromosomes[chr_num]
  #heatmap and pca
  subSNPs <- which(map$chromosome==chr  & 
                     map$physical.pos > zoom_start &
                     map$physical.pos < zoom_end
  )
  (nsnps <- length(subSNPs))
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
  
  
  ### consensus SVS colored by abundance #### 
  whichSVs <- which(SVs$chrom==chr &
                      SVs$pass2_ind #& 
                    #  (SVs$window_StartLowerCI < windows$end[i] |
                    #  SVs$window_EndUpperCI > windows$start[i])
  )
  zoom_region <- GRanges(seqnames = chr, 
                         ranges = IRanges(start = zoom_start, end = zoom_end))
  df <- SVs[whichSVs, c("chrom", "start_minus_CI_group", "end_plus_CI_group", "color")]
  localSVs <- toGRanges(df)
  kp.chrom <- plotKaryotype(plot.type=1, genome=gadmor3genome, chromosome=chr, labels.plotter=NULL,
                            zoom=zoom_region)
  kpPlotRegions(kp.chrom, data=localSVs, col=localSVs$color, border=darker("#5dade2"), r0=0, r1=1)
  kpAddBaseNumbers(kp.chrom, tick.dist=(zoom_end-zoom_start)/5, tick.len=7.5, add.units=TRUE, 
                   units="Mb", minor.tick.dist=(zoom_end-zoom_start)/10, minor.tick.len=5, cex=1)
  
  
  
    ### all pass individuals SVS colored by ecotype #### 
    whichSVs_ind <- which(SVs_ind$chrom==chr &
                        SVs_ind$pass_ind #& 
                      #  (SVs$window_StartLowerCI < windows$end[i] |
                      #  SVs$window_EndUpperCI > windows$start[i])
    )
    SVs_ind_chrom <- SVs_ind[whichSVs_ind,]
    nrow(SVs_ind_chrom)
    remove_ind <- c(which(SVs_ind_chrom$start_minus_CI_group < zoom_start &
                        SVs_ind_chrom$end_plus_CI_group < zoom_start),
                which(SVs_ind_chrom$start_minus_CI_group > zoom_end &
                        SVs_ind_chrom$end_plus_CI_group > zoom_end)
    )
    length(remove_ind)
    if(length(remove_ind)==0){remove_ind=nrow(SVs_ind_chrom)+1}
    SVs_ind_chrom <- SVs_ind_chrom[-remove_ind,]
    nrow(SVs_ind_chrom)
    SVs_ind_chrom
    if(nrow(SVs_ind_chrom)>0){
        SVs_ind_chrom <- SVs_ind_chrom[order(SVs_ind_chrom$start_minus_CI_group,
                                             SVs_ind_chrom$end_plus_CI_group,
                                             SVs_ind_chrom$pop),]
        
        plot(NULL, xlim=c(zoom_start, zoom_end), 
             ylim=c(0, nrow(SVs_ind_chrom)),
             main=c(chr, "All pass individuals in pass1 SVs"), bty="l", ylab="")
        arrows(x0=SVs_ind_chrom$start_minus_CI_group, 
               x1=SVs_ind_chrom$end_plus_CI_group,
               y0=(1:nrow(SVs_ind_chrom)), 
               lwd=2, length=0.05, code=3, angle=90, 
               col=SVs_ind_chrom$colorEcotypes
               )
        
    }    
        ### individuals SVS in pass2 SVs colored by ecotype #### 
        whichSVs_ind <- which(SVs_ind$chrom==chr &
                                SVs_ind$pass_ind &
                                SVs_ind$pass2_ind
                              #  (SVs$window_StartLowerCI < windows$end[i] |
                              #  SVs$window_EndUpperCI > windows$start[i])
        )
        SVs_ind_chrom <- SVs_ind[whichSVs_ind,]
        remove_ind <- c(which(SVs_ind_chrom$start_minus_CI_group < zoom_start &
                                SVs_ind_chrom$end_plus_CI_group < zoom_start),
                        which(SVs_ind_chrom$start_minus_CI_group > zoom_end &
                                SVs_ind_chrom$end_plus_CI_group > zoom_end)
        )
        if(length(remove_ind)==0){remove_ind=nrow(SVs_ind_chrom)+1}
        SVs_ind_chrom <- SVs_ind_chrom[-remove_ind,]
        SVs_ind_chrom <- SVs_ind_chrom[order(SVs_ind_chrom$start_minus_CI_group,
                                             SVs_ind_chrom$end_plus_CI_group,
                                             SVs_ind_chrom$pop),]
   if(nrow(SVs_ind_chrom)>0){  
        plot(NULL, xlim=c(zoom_start, zoom_end), 
             ylim=c(0, nrow(SVs_ind_chrom)),
             main=c(chr, "All pass individuals in pass1 SVs"), bty="l", ylab="")
        arrows(x0=SVs_ind_chrom$start_minus_CI_group, 
               x1=SVs_ind_chrom$end_plus_CI_group,
               y0=(1:nrow(SVs_ind_chrom)), 
               lwd=2, length=0.05, code=3, angle=90, 
               col=SVs_ind_chrom$colorEcotypes
        )
    }# end if
  dev.off()
} # end myplot

chromosomes <- levels(as.factor(map$chromosome))
chromosomes



## Loop through PCA windows and make plots #### 
## no need to rerun
 folder1 = "../results-heatmapsPCAsByMDSOutlier-100windows/"
for (j in 1:nrow(windows)){
  print(j)
    my_plots(chr_num = windows$chr_num[j], 
             zoom_start=windows$start[j], 
             zoom_end=windows$end[j], 
             folder=folder1, 
             windowID=windows$PCAwindowID[j])
  }


# Loop through chromosome level my_plots ### To do
head(gadmor3_df)
for (i in 1:nrow(gadmor3_df)){
  print(i)
  my_plots(chr_num = i, 
           zoom_start=0, 
           zoom_end=gadmor3_df$end[i]+1e06, 
           folder="../results-heatmapPCAs-chrom/", 
           windowID="")
}
 
## Chrom 1 plots ####
 my_plots(1, zoom_start <- 10*10^6, zoom_end <- 30*10^6,
          folder="../results-heatmap-otherRanges/")
  my_plots(1, zoom_start <- 12*10^6, zoom_end <- 14*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(1, zoom_start <- 15*10^6, zoom_end <- 22*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(1, zoom_start <- 15*10^6, zoom_end <- 16*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(1, zoom_start <- 17*10^6, zoom_end <- 19*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(1, zoom_start <- 20*10^6, zoom_end <- 22*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(1, zoom_start <- 25*10^6, zoom_end <- 28*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(1, zoom_start <- 28*10^6, zoom_end <- 30*10^6,
          folder="../results-heatmap-otherRanges/")
 
## Chrom 7 plots ####
 my_plots(7, zoom_start <- 15*10^6, zoom_end <- 30*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(7, zoom_start <- 16*10^6, zoom_end <- 18*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(7, zoom_start <- 23.5*10^6, zoom_end <- 24.5*10^6,
          folder="../results-heatmap-otherRanges/")
 my_plots(7, zoom_start <- 26*10^6, zoom_end <- 27*10^6,
          folder="../results-heatmap-otherRanges/")
 
 
## Chrom 12 plots ####
 my_plots(12, zoom_start <- 0*10^6, zoom_end <- 15*10^6,
          folder="../results-heatmap-otherRanges/")
 
 
 
  my_plots(chr = "NC_044051.1", zoom_start <- 17*10^6, zoom_end <- 19*10^6,
          folder="../results-heatmap-otherRanges/")
 