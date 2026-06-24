# Cod local PCA

# This shows why an enclosing circle doesnt work when run by chromosome
# I also ran it by genome and the enclosing circle doesn't work...
  # By genome the MDS2 and MDS3 are the inversions. but the grouping are weird

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



# This file would not load
#my_snps <- snp_attach("merged.f.99ind.MAF05-infos-impute.rds")

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

#seascape_full_original is the full SNP matrix - I’ve attached the structure to make it easier 
# for everyone to be on the same page. SNPs are in columns, individuals are in rows. 
# SNPs are in 0, 1, 2 format. My full script for the eLD analysis is on github:
for (chr in chromosomes) {
  chr_data_t <- t(G[,map$chromosome == chr])
  # dim(chr_data_t) SNPS in rows and individuals in columns
  
  eigenstuff <- eigen_windows(chr_data_t, win = 100, step = 10, k = 3)
  # str(eigenstuff)
  # calucate the eigenvalues and eigenvectors
  
  windist <- pc_dist(eigenstuff, npc = 3)
  # calcualte the pairwise distance between windows. npc cannot be greater than k, I think
  # str(windist)
  
  fit10d <- cmdscale(windist, eig = TRUE, k = 5)
  # get 5 MDS windows
  # str(fit10d)
  
  # The MDS coordinates can be used to find regions with "outlier" structure
  mds.coords <- fit10d$points
  pairs( mds.coords)
    # outliers for MDS 1-2
    mincirc <- getMinCircle( mds.coords[,1:2] )
    mds.corners <- corners( mds.coords[,1:2], prop=.05 )
    corner.cols <- c("red","blue","purple")
    ccols <- rep("black",nrow(mds.coords))
    for (k in 1:ncol(mds.corners)) {
      ccols[ mds.corners[,k] ] <- corner.cols[k]
    }
    plot( mds.coords[,1:2], pch=20, col=adjustcolor(ccols,0.75), 
          xlab="MDS coordinate 1", ylab="MDS coordinate 2", 
          xlim=mincirc$ctr[1]+c(-1,1)*mincirc$rad,
          ylim=mincirc$ctr[2]+c(-1,1)*mincirc$rad )
      
  
  plot(fit10d$points[,1])
  plot(fit10d$points[,2])
  plot(fit10d$points[,3])
  plot(fit10d$points[,4])
  plot(fit10d$points[,5])
  

  local_pca_results[[chr]] <- fit10d
}


## MADS CODE
# run lostruct on each chromosome with 100 SNP windows, step 20, k=2
chromosomes <- 1:10
for (chr in chromosomes) {
  chr_snps <- mut_matrix$Affx.ID[mut_matrix$Chromosome == chr]
  chr_cols <- which(colnames(seascape_full_original) %in% chr_snps)
  chr_data <- seascape_full_original[, chr_cols]
  chr_data_t <- t(chr_data)
  eigenstuff <- eigen_windows(chr_data_t, win = 100, step = 20, k = 2)
  windist <- pc_dist(eigenstuff, npc = 2)
  fit2d <- cmdscale(windist, eig = TRUE, k = 2)
  fit1d <- cmdscale(windist, eig = TRUE, k = 1)
  # max position in Mbp
  max_pos <- max(mut_matrix$Position[mut_matrix$Chromosome == chr], na.rm = TRUE) / 1e6
  window_df <- data.frame(
    Window = 1:nrow(fit2d$points),
    MDS1 = fit2d$points[, 1],
    MDS2 = fit2d$points[, 2],
    MDS1D = fit1d$points[, 1],
    Position_Mbp = seq(0, max_pos, length.out = nrow(fit2d$points)))
  local_pca_results[[chr]] <- list(
    window_df = window_df,
    windist = windist,
    fit2d = fit2d,
    fit1d = fit1d)}