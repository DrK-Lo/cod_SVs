# Cod local PCA

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

setwd("~/Desktop")


my_snps <- snp_attach("merged.f.99ind.MAF05.rds")
head(my_snps)
G <- my_snps$genotypes
dim(G)
str(G)
head(G[,1:5])
map <- my_snps$map
pops <- my_snps$fam
rm(my_snps)


# run lostruct on each chromosome with 100 SNP windows, step 20, k=2
head(map$chromosome)
chromosomes <- levels(as.factor(map$chromosome))
chromosomes
local_pca_results <- list()

# Chromsome 1 start and end based on PCA
# start at location 11299038 end at location 28292263

#seascape_full_original is the full SNP matrix - I’ve attached the structure to make it easier 
# for everyone to be on the same page. SNPs are in columns, individuals are in rows. 
# SNPs are in 0, 1, 2 format. My full script for the eLD analysis is on github:
  #chr = chromosomes[1]
  chr_data_t <- t(G[,])
  # dim(chr_data_t) SNPS in rows and individuals in columns
  
  window_size=1000
  eigenstuff <- eigen_windows(chr_data_t, win = window_size, step = 100, k = 3)
  # runs it on 1000 SNP windows. I noticed failure with 100 windows at a step of 20.
  # to do - understand how to choose k and what step means
  # str(eigenstuff)
  # calucate the eigenvalues and eigenvectors
  
  windist <- pc_dist(eigenstuff, npc = 3)
  # calcualte the pairwise distance between windows. npc cannot be greater than k, I think
  # this takes a few minutes to run on the whole genome
  # to do - understand how to choose pc_dist
  # str(windist)
  
  fit10d <- cmdscale(windist, eig = TRUE, k = 5)
  # get 5 MDS windows
  # to do - understand how to choose k
  # str(fit10d)
  
  mds.coords <- fit10d$points
  pairs( mds.coords)
  # outliers for MDS 1-2
  mincirc <- getMinCircle( mds.coords[,2:3] )
  mds.corners <- corners( mds.coords[,2:3], prop=.05 )
  corner.cols <- c("red","blue","purple")
  ccols <- rep("black",nrow(mds.coords))
  for (k in 1:ncol(mds.corners)) {
    ccols[ mds.corners[,k] ] <- corner.cols[k]
  }
  plot( mds.coords[,2:3], pch=20, col=adjustcolor(ccols,0.75), 
        xlab="MDS coordinate 1", ylab="MDS coordinate 2", 
        xlim=mincirc$ctr[1]+c(-1,1)*mincirc$rad,
        ylim=mincirc$ctr[2]+c(-1,1)*mincirc$rad )
  
  par(mfrow=c(5,1), mar=c(1,4,0,0))
  x_MB <- round(chr_map$physical.pos[seq(1,length(chr_map$physical.pos), by = window_size)]/1000000,2)
  # get position in MB for plotting
  # the window size has to match the `eigenstuff` line
  plot(fit10d$points[,1], xlab="Position", ylab="MDS1")  
  plot(fit10d$points[,2], xlab="Position", ylab="MDS2")
  plot(fit10d$points[,3], xlab="Position", ylab="MDS3")
  plot(fit10d$points[,4], xlab="Position", ylab="MDS4")
  plot(fit10d$points[,5], xlab="Position", ylab="MDS5")
  