## Consolidate V2
## KE Lotterhos
## July 2026
library(dplyr)
library(tidyr)

#setwd("/Users/k.lotterhos/Documents/GitHub/cod_SVs/consolidatev2/src")
inv <- read.csv("../../outputs/INV-dataset.csv")
head(inv)
dim(inv)

(chomosomes <- levels(as.factor(inv$CHROM)))
(PROGRAM <- levels(as.factor(inv$PROGRAM)))

head(inv[which(inv$PROGRAM=="PCA"),])

#Removing PCA regions for now
inv <- inv[-which(inv$PROGRAM=="PCA"),]
(PROGRAM <- levels(as.factor(inv$PROGRAM)))

# Remove other LGs for now
inv <- inv[which(inv$CHROM %in% chomosomes[2:24]),]
(chromosomes <- levels(as.factor(inv$CHROM)))
dim(inv)

inv <- inv %>%
  mutate(PROGRAM_col = case_when(
    PROGRAM== "DELLY" ~ "orange",
    PROGRAM == "GRIDSS" ~ "blue",
    PROGRAM == "GRIDSS (DELLY)" ~ "green",
  ))

head(inv)

## REmove 0/0 genotypes ####
table(inv$GT)
length(which(inv$GT=="0/0"))
inv <- inv[-which(inv$GT=="0/0"),]
dim(inv)

## Sanity checks ####
sum(inv$POS>inv$END, na.rm=TRUE) # should be 0 if all of the start points are less than end points
sum(is.na(inv$POS_MINUS_CI))
sum(is.na(inv$END_PLUS_CI)) # no NAs

inv$maxSize <- inv$END_PLUS_CI-inv$POS_MINUS_CI
summary(inv$maxSize)

### REvmove inversions less than 50 bp ####
hist(log10(inv$maxSize), breaks=seq(0,8,by=0.1))
sum(inv$maxSize<100, na.rm=TRUE)
inv = inv[-which(inv$maxSize<100),]
hist(log10(inv$maxSize), breaks=seq(0,8,by=0.1))

# study CI ####
hist(inv$POS_PLUS_CI-inv$POS_MINUS_CI)
hist(inv$END_PLUS_CI-inv$END_MINUS_CI)
plot(inv$POS_PLUS_CI-inv$POS_MINUS_CI, inv$END_PLUS_CI-inv$END_MINUS_CI)
plot(inv$maxSize, inv$POS_PLUS_CI-inv$POS_MINUS_CI)



## Sort inv by chromsome, pos, end ####
myorder <- order(inv$CHROM, inv$POS_MINUS_CI, inv$END_PLUS_CI, inv$PROGRAM)
inv <- inv[myorder,]
head(inv)
inv$ConsensusID <- NA


## Inspect FILTER ####
levels(as.factor(inv$FILTER))
tapply(inv$FILTER, inv$FILTER, length)

### Unique calls ####
inv$cat_unique <- paste(inv$CHROM,inv$POS_MINUS_CI,
                         inv$POS_PLUS_CI, 
                        inv$END_MINUS_CI,
                        inv$END_PLUS_CI,
                         sep="-")
head(inv)


### Binv### Being looping here, reset vars
  consensus_df = NULL
  ConsensusID_counter <- 1
  overlap <- 0.02
  inv$ConsensusID = NA
  inv_new_df = NULL
  
  hist((inv$maxSize*overlap))
  inv$window <- round(pmin(inv$maxSize*overlap,200000))
  hist(inv$window)
  
  # Extend CI
  inv$POS_PLUS_CI_window <- pmax(inv$POS_PLUS_CI,inv$POS_PLUS_CI+inv$window)
  inv$POS_MINUS_CI_window <- pmin(inv$POS_MINUS_CI,inv$POS_MINUS_CI-inv$window)
  inv$END_MINUS_CI_window <- pmin(inv$END_MINUS_CI,inv$END_MINUS_CI-inv$window)
  inv$END_PLUS_CI_window <- pmax(inv$END_PLUS_CI,inv$END_PLUS_CI+inv$window)
  
  head(inv)

### Loop through chromosomes
for (i in 1:length(chromosomes)){
  
  inv_indices <- inv$CHROM==chromosomes[i]
  print(chromosomes[i])
  inv_chrom <- inv[inv_indices,]
  print(dim(inv_chrom))
  
  pdf(paste0("../figures/",chromosomes[i],"AllCalls.pdf"), height=20, width=10)
  
  plot(NULL, xlim=c(0, max(inv_chrom$END_PLUS_CI)), ylim=c(0, nrow(inv_chrom)),
       ,main=chromosomes[i])
  arrows(x0=(inv_chrom$POS_MINUS_CI), 
         x1=inv_chrom$END_PLUS_CI,
       y0=(1:nrow(inv_chrom)), col=(inv_chrom$PROGRAM_col), xlim=c(0,30*10^6),
       main=chromosomes[i], lwd=0.1, angle=100, length=0.1, code=3)
  
  plot(NULL, xlim=c(0, max(inv_chrom$END_PLUS_CI)), ylim=c(0, nrow(inv_chrom)),
       main=chromosomes[i])
  arrows(x0=(inv_chrom$POS_MINUS_CI), 
         x1=inv_chrom$END_PLUS_CI,
           y0=(1:nrow(inv_chrom)), col=(inv_chrom$FILTER=="PASS")+1, xlim=c(0,30*10^6),
         main=chromosomes[i], lwd=0.1, angle=100, length=0.1, code=3)
  dev.off()
  
  inspect <- data.frame(table(inv_chrom$cat_unique))
  dim(inspect)
  filter <- which(inspect$Freq>5)
  length(filter)
  inspect2 <- inspect[filter,]
  #which(inv_chrom$cat_unique=="NC_044048.1-13350362-13350366-24859951-24859955")
  
  
  # visual inspection of large Chromosome 1 shows breakpoints around 
  # 13.2-13.4 MB, which is 144,000 bp (between 1-2% size of the inversion)
  # and around 24.8-24.9 or 60,000 bp
  
  while(sum(is.na(inv_chrom$ConsensusID)) > 0){
    
    (index = min(which(is.na(inv_chrom$ConsensusID))))
    
    #print(c(chromosomes[i],index,"of", nrow(inv_chrom))) #uncomment for diagnostic
    
    (focal_start_minus_CI <- inv_chrom$POS_MINUS_CI_window[index])
    (focal_start_plus_CI <- inv_chrom$POS_PLUS_CI_window[index])
    (focal_end_plus_CI <- inv_chrom$END_PLUS_CI_window[index])
    (focal_end_minus_CI <- inv_chrom$END_MINUS_CI_window[index])
    (focal_size <- inv_chrom$maxSize[index])
    
    #window = 1000
    
    ## The first group is based on a focal chromosome
    group_1 <- which(inv_chrom$POS_MINUS_CI_window<=focal_start_plus_CI &
                    (inv_chrom$POS_PLUS_CI_window)>=focal_start_minus_CI &
                    (inv_chrom$END_MINUS_CI_window)<=focal_end_plus_CI &
                    (inv_chrom$END_PLUS_CI_window)>=focal_end_minus_CI 
    )
    
    #inv_chrom[group_1,]
    dim(inv_chrom[group_1,])
    
    # The second group is based on the median breakpoint location of the first group
    (focal_start_minus_CI2 <- min(inv_chrom$POS_MINUS_CI_window[group_1], na.rm=TRUE))
    (focal_start_plus_CI2 <- max(inv_chrom$POS_PLUS_CI_window[group_1], na.rm=TRUE))
    (focal_end_plus_CI2 <- max(inv_chrom$END_PLUS_CI_window[group_1], na.rm=TRUE))
    (focal_end_minus_CI2 <- min(inv_chrom$END_MINUS_CI_window[group_1], na.rm=TRUE))
    (focal_size2 <- focal_end_plus_CI2 -focal_start_minus_CI2)
    
    group_2 <- which(inv_chrom$POS_MINUS_CI_window<=focal_start_plus_CI2 &
                       (inv_chrom$POS_PLUS_CI_window)>=focal_start_minus_CI2 &
                       (inv_chrom$END_MINUS_CI_window)<=focal_end_plus_CI2 &
                       (inv_chrom$END_PLUS_CI_window)>=focal_end_minus_CI2 
    )
 
    # Calc summary stats
    (n_calls_group <- length(group_2))
    (n_programs_group <- length(levels(as.factor(inv_chrom$PROGRAM[group_2]))))
    (n_PASS_group <- sum(inv_chrom$FILTER[group_2]=="PASS"))
    (n_ind_group <- length(levels(as.factor(inv_chrom$INDIVIDUAL[group_2]))))
    
    (start_plus_CI_group = max(inv_chrom$POS_PLUS_CI_window[group_2], na.rm=TRUE))
    (start_minus_CI_group = max(0,min(inv_chrom$POS_MINUS_CI_window[group_2]), na.rm=TRUE))
    (end_plus_CI_group = max(inv_chrom$END_PLUS_CI_window[group_2], na.rm=TRUE))
    (end_minus_CI_group = min(inv_chrom$END_MINUS_CI_window[group_2], na.rm=TRUE))
    
   
    # Save summary stats
    consensus_df_group <- data.frame(chrom = chromosomes[i],
                                      ConsensusID=ConsensusID_counter,
                                     n_calls_group,
                                     n_programs_group,
                                     n_PASS_group,
                                     n_ind_group,
                                     start_plus_CI_group,
                                     start_minus_CI_group,
                                     end_plus_CI_group,
                                     end_minus_CI_group)
    consensus_df <- rbind(consensus_df,consensus_df_group)
    
    # Give them the same ID
    dim(inv_chrom[group_2,])
    inv_chrom$ConsensusID[group_2] <- ConsensusID_counter
    
    ConsensusID_counter <- ConsensusID_counter + 1
  } #end while loop
  
  inv_new_df <- rbind(inv_new_df,inv_chrom)
    
} # end chrom loop ####


## AFter looping - Summary stats for consensus calls and PASS criteria ####

head(consensus_df)
dim(consensus_df)
  plot(consensus_df$n_programs_group, consensus_df$n_ind_group)

  consensus_df$pass1 <- consensus_df$n_programs_group >= 3 &
                      consensus_df$n_PASS_group >= 1 &
                      consensus_df$n_ind_group >= 15

  table(consensus_df$chrom)
  table(consensus_df$chrom,consensus_df$pass1)
  sum(consensus_df$pass1)


### Filtered individuals DF ####
dim(inv_new_df)
ind_filtered <- merge(inv_new_df[,c("INDIVIDUAL","PROGRAM","GT","FILTER","ConsensusID")], consensus_df, all.x=TRUE)
dim(ind_filtered) # same dimensions as pre-merge
head(ind_filtered)

ind_filtered <- ind_filtered[which(ind_filtered$pass1),]
dim(ind_filtered)
head(ind_filtered)

## Count number of calls supporting that INV within the individual
ind_filtered$InvIndID <- paste(ind_filtered$INDIVIDUAL, ind_filtered$ConsensusID, sep="__")
n_calls_ind <- data.frame(table(ind_filtered$InvIndID))
colnames(n_calls_ind) <- c("InvIndID", "n_calls_ind")
head(n_calls_ind)
dim(ind_filtered)
ind_filtered <- merge(ind_filtered, n_calls_ind, all.x=TRUE)
dim(ind_filtered)

## Count number of programs supporting that INV within the individual
n_programs_ind <- table(ind_filtered$InvIndID, ind_filtered$PROGRAM)
n_programs_ind <- n_programs_ind > 0
n_programs_ind2 <- data.frame(InvIndID=rownames(n_programs_ind),
                              n_programs_ind=rowSums(n_programs_ind))
head(n_programs_ind2)
ind_filtered <- merge(ind_filtered, n_programs_ind2, all.x=TRUE)
dim(ind_filtered)
head(ind_filtered)

## Count number of PASS supporting that INV within the individual
n_pass_ind <- table(ind_filtered$InvIndID, ind_filtered$FILTER)
head(n_pass_ind)
n_pass_ind2 <- data.frame(InvIndID=rownames(n_pass_ind),
                          n_pass_ind=n_pass_ind[,"PASS"])
head(n_pass_ind2)
ind_filtered <- merge(ind_filtered, n_pass_ind2, all.x=TRUE)
dim(ind_filtered)
head(ind_filtered)


#### Order and filter individuals ####
ind_filtered <- ind_filtered[order(ind_filtered$INDIVIDUAL,
                                   ind_filtered$ConsensusID,
                                   ind_filtered$PROGRAM),]
head(ind_filtered,50)

ind_filtered$dup <- duplicated(ind_filtered$InvIndID)
sum(ind_filtered$dup)

### Final individual filtered ####
ind_filtered_final <- ind_filtered[!ind_filtered$dup,]
dim(ind_filtered_final)


ind_filtered_final$pop.names <- substr(ind_filtered_final$INDIVIDUAL,0,4)
ind_filtered_final$samp_full <- substr(ind_filtered_final$INDIVIDUAL,6,10)
tail(ind_filtered_final) 
 
ind_filtered_final <- ind_filtered_final[,-which(colnames(ind_filtered_final) %in% c("GT", "PROGRAM", "FILTER", "dup"))]
head(ind_filtered_final)

### Read in population and sample data and annotate the individual table
popdata <- read.table("../../data/1-PopsColors.txt", header=TRUE)
sampledata <- readRDS("../../data/2-Samples.rds")

popdata
dim(ind_filtered_final)
ind_filtered_final2 <- merge(ind_filtered_final, popdata[,c("pop.names","pop","colorEcotypes")], all.x=TRUE)
dim(ind_filtered_final2)
head(ind_filtered_final2)

head(sampledata)
ind_filtered_final3 <- merge(ind_filtered_final2, 
                             sampledata[,c("samp_full",
                                           "LG1cluster_nameB",
                                           "LG2cluster_nameB",
                                           "LG7cluster_nameB",
                                           "LG12cluster_nameB")], all.x=TRUE)
dim(ind_filtered_final3)
head(ind_filtered_final3)


ind_filtered_final3 <- ind_filtered_final3[order(ind_filtered_final3$INDIVIDUAL,
                                   ind_filtered_final3$ConsensusID),]



## Individual pass criteria ####
ind_filtered_final3$pass_ind <- ind_filtered_final3$n_programs_ind >= 3 & 
  ind_filtered_final3$n_pass_ind > 0


#number of inversions per individual ####
inv_per_ind <- table(ind_filtered_final3$INDIVIDUAL)
hist(inv_per_ind)
summary(as.numeric(inv_per_ind))


inv_per_ind <- table(ind_filtered_final3$INDIVIDUAL[which(ind_filtered_final3$pass_ind)])
hist(inv_per_ind)
summary(as.numeric(inv_per_ind))

str(ind_filtered_final3)

## Calculate individual level pass stats for each Consensus ID ####
  ind_pass <- table(ind_filtered_final3$ConsensusID, ind_filtered_final3$pass_ind)
  ind_pass <- data.frame(as.matrix(ind_pass))
  head(ind_pass)
  names(ind_pass) <- c("ConsensusID", "pass", "n_ind")
  head(ind_pass)
  ind_pass <- data.frame(pivot_wider(ind_pass, names_from = pass, values_from = n_ind))
  head(ind_pass)
  colnames(ind_pass) <- c("ConsensusID","n_NOpass_ind_group", "n_PASS_ind_group")
  #write.csv(ind_pass, "../data_outputs/INV-dataset-ConsensusID_unique-IndLevelPassStats.csv")
  hist(ind_pass[,"n_PASS_ind_group"])
  ind_pass$pass2_ind <- ind_pass$n_PASS_ind > 5
  head(ind_pass)
  
  dim(consensus_df)
  consensus_df2 <- merge(consensus_df, ind_pass, by="ConsensusID", all.x=TRUE)
  dim(consensus_df2)
  head(consensus_df2)

  # inspect 
  consensus_df2$maxSize <- consensus_df2$end_plus_CI_group-consensus_df2$start_minus_CI_group
  LG12 <- consensus_df2[consensus_df2$chrom==chromosomes[12],]
  LG12 <- LG12[which(LG12$end_minus_CI_group<15*10^6 & LG12$maxSize>8*10^6),]
  # LG12 is tricky. The consensus ID for the likely inversion is 11797.
  # It starts near 0 and ends near 15, and 21 individuals share it. 
  # However, within an individual, only one individual passes.
  # There are others, but no passing
  
## Loop through chrom and make plots ####
for (i in 1:length(chromosomes)){
  
  print(chromosomes[i])
  
  consensus_df_chrom <- consensus_df2[consensus_df$chrom==chromosomes[i],]
  consensus_df_chrom <- consensus_df_chrom[order(consensus_df_chrom$chrom, 
                                                 consensus_df_chrom$start_minus_CI_group,
                                                 consensus_df_chrom$end_plus_CI_group),]

  consensus_pass1 <- consensus_df_chrom[which(consensus_df_chrom$pass1),]
  consensus_pass2 <- consensus_df_chrom[which(consensus_df_chrom$pass2_ind),]

pdf(paste0("../figures/",chromosomes[i],"ConsensusCalls.pdf"), height=20, width=10)

#hist(consensus_df_chrom$start_minus_CI_group/10^6, breaks=0:60)
#hist(consensus_df_chrom$end_plus_CI_group/10^6, breaks=0:60)

plot(NULL, xlim=c(0, max(consensus_df_chrom$end_plus_CI_group)), 
     ylim=c(0, nrow(consensus_df_chrom)),
     main=chromosomes[i])
arrows(x0=(consensus_df_chrom$start_plus_CI_group+
             consensus_df_chrom$start_minus_CI_group)/2, 
       x1=(consensus_df_chrom$end_plus_CI_group+
             consensus_df_chrom$end_minus_CI_group)/2,
       y0=(1:nrow(consensus_df_chrom)), 
       lwd=consensus_df_chrom$n_ind_group/20, col="red", 
       angle=100, length=0.1, code=3)


#hist(consensus_pass$start_minus_CI_group/10^6, breaks=0:60)
#hist(consensus_pass$end_plus_CI_group/10^6, breaks=0:60)

plot(NULL, xlim=c(0, max(consensus_df_chrom$end_plus_CI_group)), 
     ylim=c(0, nrow(consensus_pass1)),
     main=c(chromosomes[i], "pass1"))
arrows(x0=(consensus_pass1$start_plus_CI_group+
             consensus_pass1$start_minus_CI_group)/2, 
       x1=(consensus_pass1$end_plus_CI_group+
             consensus_pass1$end_minus_CI_group)/2,
       y0=(1:nrow(consensus_pass1)), lwd=consensus_pass1$n_ind_group/20+1, 
       col="red", angle=100, length=0.1, code=3)

plot(NULL, xlim=c(0, max(consensus_df_chrom$end_plus_CI_group)), 
     ylim=c(0, nrow(consensus_pass2)),
     main=c(chromosomes[i], "pass2"))
arrows(x0=(consensus_pass2$start_plus_CI_group+
             consensus_pass2$start_minus_CI_group)/2, 
       x1=(consensus_pass2$end_plus_CI_group+
             consensus_pass2$end_minus_CI_group)/2,
       y0=(1:nrow(consensus_pass2)), lwd=consensus_pass2$n_ind_group/20+1, 
       col="red", angle=100, length=0.1, code=3)
dev.off()

}# end loop through chromosomes

  
### Write tables ####
  write.csv(consensus_df2,"../data_outputs/INV-dataset-ConsensusID_unique.csv")
  dim(consensus_df2)
  hist(consensus_df2$ConsensusID, breaks=seq(0,30000, by=1000))
  
  dim(inv)
  dim(inv_new_df)
  sum(is.na(inv_new_df$ConsensusID)) # sanity check, should be 0
  hist(inv_new_df$ConsensusID, breaks=seq(0,30000, by=1000)) #sanity check, make sure filtered IDs span range
  write.csv(inv_new_df,gzfile("../data_outputs/INV-dataset-ConsensusID-unfiltered-individuals.csv.gz"))
  
  
## Write Individual file with filtering ####
  head(ind_filtered_final3)
  dim(ind_filtered_final3)
  head(ind_pass)
  ind_filtered_final4 <- merge(ind_filtered_final3, ind_pass, by="ConsensusID", all.x=TRUE) #merge in pass 2 criteria
  dim(ind_filtered_final4)
  head(ind_filtered_final4)


write.csv(ind_filtered_final4,"../data_outputs/INV-dataset-ConsensusID-Ind-Filtered.csv")

## Loop through chrom and make individual level plots ####
for (i in 1:length(chromosomes)){
  
  print(chromosomes[i])
  
  ind_filtered_final3_chrom <- ind_filtered_final3[ind_filtered_final3$chrom==chromosomes[i],]
  ind_filtered_final3_chrom <- ind_filtered_final3_chrom[order(ind_filtered_final3_chrom$chrom, 
                                                 ind_filtered_final3_chrom$start_minus_CI_group,
                                                 ind_filtered_final3_chrom$end_plus_CI_group,
                                                 ind_filtered_final3_chrom$pop),]
  head(ind_filtered_final3_chrom)
  ind_filtered_final3_chrom_filtered <- ind_filtered_final3_chrom[which(ind_filtered_final3_chrom$pass_ind),]
  
  dim(ind_filtered_final3_chrom )
  dim(ind_filtered_final3_chrom_filtered)
  
  pdf(paste0("../figures/",chromosomes[i],"ConsensusCalls_IndLevel.pdf"), height=20, width=10)
  
  plot(NULL, xlim=c(0, max(ind_filtered_final3_chrom$end_plus_CI_group)), 
       ylim=c(0, nrow(ind_filtered_final3_chrom)),
       main=chromosomes[i])
  arrows(x0=(ind_filtered_final3_chrom$start_plus_CI_group+
               ind_filtered_final3_chrom$start_minus_CI_group)/2, 
         x1=(ind_filtered_final3_chrom$end_plus_CI_group+
               ind_filtered_final3_chrom$end_minus_CI_group)/2,
         y0=(1:nrow(ind_filtered_final3_chrom)), 
         col=as.character(ind_filtered_final3_chrom$colorEcotypes),
         code=3, angle=100, length=0.1)
  
  plot(NULL, xlim=c(0, max(ind_filtered_final3_chrom$end_plus_CI_group)), 
       ylim=c(0, nrow(ind_filtered_final3_chrom_filtered)),
       main=paste(chromosomes[i], "filtered 3 programs; 1 pass"))
   arrows(x0=(ind_filtered_final3_chrom_filtered$start_plus_CI_group+
               ind_filtered_final3_chrom_filtered$start_minus_CI_group)/2, 
         x1=(ind_filtered_final3_chrom_filtered$end_plus_CI_group+
               ind_filtered_final3_chrom_filtered$end_minus_CI_group)/2,
         y0=(1:nrow(ind_filtered_final3_chrom_filtered)), 
         col=as.character(ind_filtered_final3_chrom_filtered$colorEcotypes),
         angle=100, length=0.1, code=3)

  dev.off()
  
}# end loop through chromosomes
   
