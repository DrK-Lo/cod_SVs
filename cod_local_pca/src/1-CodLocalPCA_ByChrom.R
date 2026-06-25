# Cod local PCA
# Author: KE Lotterhos
# Date: June 2026

#update.packages(ask = FALSE, checkBuilt = TRUE)
#install.packages("devtools")
#install.packages("rlang", dependencies=TRUE)
#install.packages("shotGroups")
library(rlang)
#devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
library(bigstatsr)
library(bigsnpr)
library(Matrix)
library(dplyr)
library(LEA)

setwd("~/Documents/GitHub/cod_SVs/cod_local_pca/src")

my_snps <- snp_attach("~/Desktop/codGenotypeData/merged.f.99ind.MAF05.rds")
# this file is too large for the GitHub

head(my_snps)

# Get genotypes and remove duplicate individuals
G <- my_snps$genotypes
G <- G[,] # make G a matrix

dupIndex <- which(duplicated(my_snps$fam$sample.ID)) # find the duplicate
my_snps$fam$sample.ID[dupIndex]
remove = which(my_snps$fam$sample.ID == my_snps$fam$sample.ID[dupIndex])[1]
remove

# Check for missing data at the individual level
missing <- rowSums(is.na(G))
summary(missing) # one individual is missing a lot!
remove <- c(remove, which(missing == max(missing)))
remove

# Remove the duplicates and missing data
G <- G[-remove,]
dim(G)
str(G)
head(G[,1:5])
map <- my_snps$map
pops <- my_snps$fam[-remove,]

# double check missing data
missing <- rowSums(is.na(G))
summary(missing) # 50K SNPs is about 5%, which is OK

rm(my_snps)



# run lostruct on each chromosome with 100 SNP windows, step 20, k=2
head(map$chromosome)
chromosomes <- levels(as.factor(map$chromosome))
chromosomes

#### prep for karyoploter ####
#Create custom genome
gadmor3_df <- data.frame(chr=c("NC_044048.1", "NC_044049.1", "NC_044050.1", "NC_044051.1", "NC_044052.1", "NC_044053.1", "NC_044054.1", "NC_044055.1", "NC_044056.1", "NC_044057.1", "NC_044058.1", "NC_044059.1", "NC_044060.1", "NC_044061.1", "NC_044062.1", "NC_044063.1", "NC_044064.1", "NC_044065.1", "NC_044066.1", "NC_044067.1", "NC_044068.1", "NC_044069.1", "NC_044070.1"), start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), end=c(30875876, 28732775, 30954429, 43798135, 25300426, 27762770, 34137969, 29710654, 26487948, 27234273, 30713045, 30948897, 28829685, 29586942, 28657694, 34794352, 21723002, 24902675, 22015597, 24843429, 22358821, 23744039, 25242006))
head(gadmor3_df)
gadmor3genome <- toGRanges(gadmor3_df)
head(gadmor3genome)

# Chromsome 1 start and end based on PCA
# start at location 11299038 end at location 28292263

#seascape_full_original is the full SNP matrix - I’ve attached the structure to make it easier 
# for everyone to be on the same page. SNPs are in columns, individuals are in rows. 
# SNPs are in 0, 1, 2 format. My full script for the eLD analysis is on github:
local_pca_results <- NULL
maxPCAwindowID <- 0
dir.create(paste0("../results-MDSbyChrom-", window_size,"windows"))
for (chr in chromosomes) {
  #chr = chromosomes[1] #uncomment for testing
  chr_num <- as.numeric(substr(chr,8,9))-47
  print(chr_num)
  chr_data_t <- t(G[,which(map$chromosome == chr)])
  chr_map <- map[which(map$chromosome == chr),]
  # dim(chr_data_t) SNPS in rows and individuals in columns
  
  window_size=100
  eigenstuff <- eigen_windows(chr_data_t, win = window_size, step = 10, k = 3)
  # It runs on 1000 SNP windows with a step of 20. I noticed failure with 100 windows at a step of 20.
  # to do - understand how to choose k and what step means
  # str(eigenstuff)
  # calucate the eigenvalues and eigenvectors
  
  windist <- pc_dist(eigenstuff, npc = 3)
  # calcualte the pairwise distance between windows. npc cannot be greater than k, I think
  # to do - understand how to choose  npc
  # str(windist)
  #heatmap(windist)
  sum(is.na(windist)) # should be non-zero
  
  fit10d <- cmdscale(windist, eig = TRUE, k = 5)
  # get 5 MDS windows
  # str(fit10d)
  
  # get position in MB for plotting
  # the window size has to match the `eigenstuff` line
  x <- chr_map$physical.pos[seq(1,length(chr_map$physical.pos),by = window_size)]
  x_MB <- round( x/1000000,2)
  

  # Outlier regions are based on those that are num_sds SD more extreme than the median
  num_sds <- 3
  outlier_MDS1 <- (fit10d$points[,1] < median(fit10d$points[,1])-num_sds*sd(fit10d$points[,1])) | 
                  (fit10d$points[,1] > median(fit10d$points[,1])+num_sds*sd(fit10d$points[,1]))
  outlier_MDS2 <- (fit10d$points[,2] < median(fit10d$points[,2])-num_sds*sd(fit10d$points[,2])) | 
                  (fit10d$points[,2] > median(fit10d$points[,2])+num_sds*sd(fit10d$points[,2]))
  outlier_MDS3 <- (fit10d$points[,3] < median(fit10d$points[,3])-num_sds*sd(fit10d$points[,3])) | 
                  (fit10d$points[,3] > median(fit10d$points[,3])+num_sds*sd(fit10d$points[,3]))
  outlier_MDS4 <- (fit10d$points[,4] < median(fit10d$points[,4])-num_sds*sd(fit10d$points[,4])) | 
                  (fit10d$points[,4] > median(fit10d$points[,4])+num_sds*sd(fit10d$points[,4]))
  outlier_MDS5 <- (fit10d$points[,5] < median(fit10d$points[,5])-num_sds*sd(fit10d$points[,5])) | 
                  (fit10d$points[,5] > median(fit10d$points[,5])+num_sds*sd(fit10d$points[,5]))
  
  
  if (chr == "NC_044048.1"){ # hard code Chrom 1 as outlier
    outlier_MDS1 <- rep(FALSE, length=nrow(fit10d$points))
    if (window_size==1000){
      outlier_MDS1[which(fit10d$points[,1] > 
                           median(fit10d$points[1:40,1])+1.96*sd(fit10d$points[1:40,1])
      )] <- TRUE
    }
    if (window_size ==100){
      outlier_MDS1[which(fit10d$points[,1] > 
                           median(fit10d$points[1:400,1])+1.96*sd(fit10d$points[1:400,1])
      )] <- TRUE
    }

  }
  
  if (chr == "NC_044049.1"){ # hard code Chrom 2 as outlier
    outlier_MDS1 <- rep(FALSE, length=nrow(fit10d$points))
    if (window_size ==100){
      outlier_MDS1[which(fit10d$points[,1] > 
                           median(fit10d$points[400:nrow(fit10d$points),1])+num_sds*
                           sd(fit10d$points[400:nrow(fit10d$points),1])
      )] <- TRUE
    }
  }
  if (chr == "NC_044054.1"){ # hard code chromosome 7 as outlier
    outlier_MDS1 <- rep(FALSE, length=nrow(fit10d$points))
    if (window_size ==1000){
      outlier_MDS1[which(fit10d$points[,1] > 
                           median(fit10d$points[c(1:50,150:169),1])+num_sds
                         *sd(fit10d$points[c(1:50,150:169),1])
      )] <- TRUE
    }
    if (window_size ==100){
      outlier_MDS1[which(fit10d$points[,1] > 
                           median(fit10d$points[c(1:500,1500:1690),1])+num_sds*
                           sd(fit10d$points[c(1:50,1500:1690),1])
      )] <- TRUE
    }
  }
  if (chr == "NC_044059.1"){ # hard code chromosome 12 as outliers
    outlier_MDS1 <- rep(FALSE, length=nrow(fit10d$points))
    outlier_MDS1[which(fit10d$points[,1] > 
                         median(fit10d$points[700:1230,1])+num_sds*
                         sd(fit10d$points[700:1230,1])
    )] <- TRUE
  }
  
  # outliers based on quantiles
 # lower_quantile=0.025
  # upper_quantile=0.975
    # I found this didn't work well because the number of outliers is fixed based on the number of windows
  # outlier_MDS1 <- (fit10d$points[,1] < quantile(fit10d$points[,1],lower_quantile)) | (fit10d$points[,1] > quantile(fit10d$points[,1],upper_quantile))
  # outlier_MDS2 <- (fit10d$points[,2] < quantile(fit10d$points[,2],lower_quantile)) | (fit10d$points[,2] > quantile(fit10d$points[,2],upper_quantile))
  # outlier_MDS3 <- (fit10d$points[,3] < quantile(fit10d$points[,3],lower_quantile)) | (fit10d$points[,3] > quantile(fit10d$points[,3],upper_quantile))
  # outlier_MDS4 <- (fit10d$points[,4] < quantile(fit10d$points[,4],lower_quantile)) | (fit10d$points[,4] > quantile(fit10d$points[,4],upper_quantile))
  # outlier_MDS5 <- (fit10d$points[,5] < quantile(fit10d$points[,5],lower_quantile)) | (fit10d$points[,5] > quantile(fit10d$points[,5],upper_quantile))
  # 

  
    PCA_outlier <- outlier_MDS1 | outlier_MDS2 | outlier_MDS3 | outlier_MDS4 | outlier_MDS5
    
    PCAwindowID_thischrom <- consecutive_id(PCA_outlier)*as.numeric(PCA_outlier)
      # this gives an ID to non-outliers
    
    PCAwindowID_thischrom[which(PCAwindowID_thischrom==0)] <- NA
      # give the non-outliers an NA
    
    if (min(PCAwindowID_thischrom, na.rm=TRUE)==1){PCAwindowID_thischrom <- (PCAwindowID_thischrom+1)/2}  
    if (min(PCAwindowID_thischrom, na.rm=TRUE)==2){PCAwindowID_thischrom <- (PCAwindowID_thischrom)/2}
    
      # make the IDs continuous
    
    PCAwindowID <- PCAwindowID_thischrom + maxPCAwindowID
    maxPCAwindowID <- max(PCAwindowID, na.rm=TRUE)
    
    # create a vector with one label per window for plotting
    thisChromRow=1:length(outlier_MDS1)
    floor <- floor(tapply(thisChromRow,PCAwindowID,mean))
    PCAwindowID_plot <- rep(NA, length(PCAwindowID))
    PCAwindowID_plot[floor] <- PCAwindowID[floor]
    
  this_chrom <- data.frame(chromosome=chr, chr_num,
             start_pos=x[1:nrow(fit10d$points)], 
             end_pos=x[2:length(x)]-1, 
             start_MB = x_MB[1:nrow(fit10d$points)],
             outlier_MDS1,
             outlier_MDS2,
             outlier_MDS3,
             outlier_MDS4,
             outlier_MDS5,
             PCA_outlier,
             PCAwindowID,
             thisChromRow,
             PCAwindowID_plot
             )
  head(this_chrom)
  local_pca_results <- bind_rows(local_pca_results,this_chrom) # add each chromosome to the data frame

    
  #png(paste0("results-MDSbyChrom/",chr,"_localPCA_",chr_num,"_",window_size,"SNPwindows.png"),res=500, width=6, units="in", height=7)
  pdf(paste0("../results-MDSbyChrom-",window_size,"windows/",chr,"_localPCA_",chr_num,"_",window_size,"SNPwindows.pdf"),width=6, height=7)
  par(mfrow=c(6,1), mar=c(1,4,0.4,0), oma=c(4,0,4,0))
  PCA_outlier2 <- PCA_outlier
  PCA_outlier2[which(PCA_outlier==FALSE)]=NA

  plot(this_chrom$start_MB,rep(1, times=nrow(fit10d$points)),  ylim=c(0.9,1.1), col="grey", ylab="",xaxt="n", yaxt="n", bty="n", pch=19)
  points(this_chrom$start_MB, PCA_outlier2, lwd=5,  col="red")
    text(this_chrom$start_MB[!is.na(this_chrom$PCAwindowID_plot)], y=rep(c(1.05, 1.08), times=nrow(fit10d$points)/2), 
         labels = this_chrom$PCAwindowID_plot[!is.na(this_chrom$PCAwindowID_plot)], cex=0.5, srt=90)
 
    
  plot(this_chrom$start_MB, fit10d$points[,1], xlab="Position", ylab="MDS1", col=outlier_MDS1+1, pch=outlier_MDS1+20, xaxt="n", bty="l", las=1, cex=1.5)  
  plot(this_chrom$start_MB, fit10d$points[,2], xlab="Position", ylab="MDS2", col=outlier_MDS2+1, pch=outlier_MDS2+20, xaxt="n", bty="l", las=1, cex=1.5) 
  plot(this_chrom$start_MB, fit10d$points[,3], xlab="Position", ylab="MDS3", col=outlier_MDS3+1, pch=outlier_MDS3+20, xaxt="n", bty="l", las=1, cex=1.5) 
  plot(this_chrom$start_MB, fit10d$points[,4], xlab="Position", ylab="MDS4", col=outlier_MDS4+1, pch=outlier_MDS4+20, xaxt="n", bty="l", las=1, cex=1.5)
  plot(this_chrom$start_MB, fit10d$points[,5], ylab="MDS5", col=outlier_MDS5+1, pch=outlier_MDS5+20,  bty="l", las=1, cex=1.5, xlab="Position (MB)")
  mtext(paste0("Chromosome ", chr_num," (" ,chr, ")"), side=3, outer=TRUE, adj=0.1)
  mtext("Position (MB)", side=1, outer=TRUE, line=2)

  dev.off()
  #dev.off()
  
  pdf(paste0("../results-MDSbyChrom-",window_size,"windows/",chr,"_localPCA_",chr_num,"_",window_size,"SNPwindows_KaryoPlot.pdf"),width=6, height=2)
  start <- tapply(this_chrom$start_pos,this_chrom$PCAwindowID, min)
  end <- tapply(this_chrom$end_pos,this_chrom$PCAwindowID, max)
  ranges <- toGRanges(data.frame(Chrom=chr, start, end))
  par(mfrow=c(1,1), mar=c(1,4,0.4,0), oma=c(0,0,0,0))
  kp.chrom <- plotKaryotype(plot.type=1, genome=gadmor3genome, chromosome=chr, labels.plotter=NULL)
  kpPlotRegions(kp.chrom, data=ranges, col="#5dade2", border=darker("#5dade2"), r0=0, r1=1)
  kpAddBaseNumbers(kp.chrom, tick.dist=10000000, tick.len=7.5, add.units=TRUE, units="Mb", minor.tick.dist=5000000, minor.tick.len=5, cex=1)
  seqlevels(kp.chrom$genome)
  kpText(kp.chrom, chr=chr, data=ranges, x = (ranges$end-ranges$start)/2,
         y=c(1.1,1.3), labels = names(start), srt=90, cex=0.4)
  dev.off()
}  



write.table(local_pca_results, 
            paste0("../data_outputs/local_pca_results-",window_size,"SNPwindow.txt"))

# Mean window size in bp
head(local_pca_results)
dim(local_pca_results)
summary(local_pca_results$end_pos-local_pca_results$start_pos)
  # 100 base pair windows give mean window size of 23568 and median of 19912 bases
dim(G)

### Get a list of unique PCA windows and their start and end points
start <- tapply(local_pca_results$start_pos,local_pca_results$PCAwindowID, min)
end <- tapply(local_pca_results$end_pos,local_pca_results$PCAwindowID, max)
numWin <- tapply(local_pca_results$PCAwindowID,local_pca_results$PCAwindowID, length)
AllPCAwindows <- data.frame(PCAwindowID=as.numeric(names(start)), start, end, numSNPs = numWin*window_size)
ChomIDS <- local_pca_results[,c("chromosome","PCAwindowID","chr_num" )]
ChomIDS <- ChomIDS %>% distinct(PCAwindowID, .keep_all = TRUE)
str(ChomIDS)

summary(AllPCAwindows$end - AllPCAwindows$start)
summary(AllPCAwindows$numSNPs)
  
str(AllPCAwindows)
str(local_pca_results)

AllPCAwindows2 <- left_join(AllPCAwindows, ChomIDS)
write.table(AllPCAwindows2, paste0("../data_outputs/local_pca_results_UniquePCAwindows-",
                                   window_size,"SNPwindow.txt"))
