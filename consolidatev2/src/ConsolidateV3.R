## Consolidate V2
## KE Lotterhos
## July 2026
library(dplyr)

setwd("/Users/k.lotterhos/Documents/GitHub/cod_SVs/consolidatev2/src")
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
    PROGRAM== "DELLY" ~ "red",
    PROGRAM == "GRIDSS" ~ "blue",
    PROGRAM == "GRIDSS (DELLY)" ~ "purple",
  ))

head(inv)

## Sanity checks
sum(inv$POS>inv$END, na.rm=TRUE) # should be 0 if all of the start points are less than end points
sum(is.na(inv$POS_MINUS_CI))
sum(is.na(inv$END_PLUS_CI))

inv$maxSize <- inv$END_PLUS_CI-inv$POS_MINUS_CI
summary(inv$maxSize)
inv <- inv[-which(inv$maxSize<100),]
summary(inv$maxSize)
hist(log10(inv$maxSize))

# study CI ####
hist(inv$POS_PLUS_CI-inv$POS_MINUS_CI)
hist(inv$END_PLUS_CI-inv$END_MINUS_CI)
plot(inv$POS_PLUS_CI-inv$POS_MINUS_CI, inv$END_PLUS_CI-inv$END_MINUS_CI)
plot(inv$maxSize, inv$POS_PLUS_CI-inv$POS_MINUS_CI)



## Sort inv by chromsome,POS_PLUS_CI## Sort inv by chromsome, pos, end
myorder <- order(inv$CHROM, inv$POS_MINUS_CI, inv$END_PLUS_CI, inv$PROGRAM)
inv <- inv[myorder,]
head(inv)
inv$ConsensusID <- NA


## Inspect FILTER
levels(as.factor(inv$FILTER))
tapply(inv$FILTER, inv$FILTER, length)

### Unique inversions
inv$cat_unique <- paste(inv$CHROM,inv$POS_MINUS_CI,
                         inv$POS_PLUS_CI, 
                        inv$END_MINUS_CI,
                        inv$END_PLUS_CI,
                         sep="-")
head(inv)



## Extend CI


### Binv### Being looping here, reset vars
consensus_df = NULL
ConsensusID_counter <- 1
overlap <- 0.02
inv$ConsensusID = NA
inv_new_df = NULL

hist((inv$maxSize*overlap))
inv$window <- round(pmin(inv$maxSize*overlap,200000))
hist(inv$window)
inv$POS_PLUS_CI_window <- pmax(inv$POS_PLUS_CI,inv$POS_PLUS_CI+inv$window)
inv$POS_MINUS_CI_window <- pmin(inv$POS_MINUS_CI,inv$POS_MINUS_CI-inv$window)
inv$END_MINUS_CI_window <- pmin(inv$END_MINUS_CI,inv$END_MINUS_CI-inv$window)
inv$END_PLUS_CI_window <- pmax(inv$END_PLUS_CI,inv$END_PLUS_CI+inv$window)

head(inv)

### Loop through chromosomes
for (i in 1:length(chromosomes)){
  
  inv_indices <- inv$CHROM==chromosomes[i]
  inv_chrom <- inv[inv_indices,]
  print(dim(inv_chrom))
  
  pdf(paste0("../figures/",chromosomes[i],"AllCalls.pdf"), height=20, width=10)
  
  plot(NULL, xlim=c(0, max(inv_chrom$END_PLUS_CI)), ylim=c(0, nrow(inv_chrom)),
       ,main=chromosomes[i])
  arrows(x0=(inv_chrom$POS_MINUS_CI), 
         x1=inv_chrom$END_PLUS_CI,
       y0=(1:nrow(inv_chrom)), col=(inv_chrom$PROGRAM_col), xlim=c(0,30*10^6),
       main=chromosomes[i], lwd=0.1)
  
  plot(NULL, xlim=c(0, max(inv_chrom$END_PLUS_CI)), ylim=c(0, nrow(inv_chrom)),
       main=chromosomes[i])
  arrows(x0=(inv_chrom$POS_MINUS_CI), 
         x1=inv_chrom$END_PLUS_CI,
           y0=(1:nrow(inv_chrom)), col=(inv_chrom$FILTER=="PASS")+1, xlim=c(0,30*10^6),
         main=chromosomes[i], lwd=0.1)
  dev.off()
  
  inspect <- data.frame(table(inv_chrom$cat_unique))
  dim(inspect)
  filter <- which(inspect$Freq>15)
  length(filter)
  inspect[filter,]
  #which(inv_chrom$cat_unique=="NC_044048.1-13350362-13350366-24859951-24859955")
  
  
  # visual inspection of large Chromosome 1 shows breakpoints around 
  # 13.2-13.4 MB, which is 144,000 bp (between 1-2% size of the inversion)
  # and around 24.8-24.9 or 60,000 bp
  
  while(sum(is.na(inv_chrom$ConsensusID)) > 0){
    
    (index = min(which(is.na(inv_chrom$ConsensusID))))
    
    print(c(index,"of", nrow(inv_chrom)))
    
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
                                     end_minus_CI_group,
                                      pass=NA)
    consensus_df <- rbind(consensus_df,consensus_df_group)
    
    # Give them the same ID
    dim(inv_chrom[group_2,])
    inv_chrom$ConsensusID[group_2] <- ConsensusID_counter
    
    ConsensusID_counter <- ConsensusID_counter + 1
  } #end while loop
  
  inv_new_df <- rbind(inv_new_df,inv_chrom)
    
} # end chrom loop ####


## Summary stats for consensus calls 

head(consensus_df)
dim(consensus_df)

  consensus_df$pass <- consensus_df$n_programs_group >= 2 &
                      consensus_df$n_PASS_group >= 1 &
                      consensus_df$n_ind_group > 15

  table(consensus_df$chrom)
  table(consensus_df$chrom,consensus_df$pass)

## Loop through chrom and make plots
for (i in 1:length(chromosomes)){
  
  print(chromosomes[i])
    
  consensus_df_chrom <- consensus_df[consensus_df$chrom==chromosomes[i],]
  consensus_df_chrom <- consensus_df_chrom[order(consensus_df_chrom$chrom, 
                                         consensus_df_chrom$start_minus_CI_group,
                                         consensus_df_chrom$end_plus_CI_group),]
  mypass <- which(consensus_df_chrom$pass)
  npass <- length(mypass)

  consensus_pass <- consensus_df_chrom[mypass,]
  dim(consensus_pass)
  consensus_pass <- consensus_pass[order(consensus_pass$chrom, 
                                     consensus_pass$start_minus_CI_group,
                                     consensus_pass$end_plus_CI_group),]


pdf(paste0("../figures/",chromosomes[i],"ConsensusCalls.pdf"), height=20, width=10)
  
 hist(consensus_df_chrom$start_minus_CI_group/10^6, breaks=0:60)
 hist(consensus_df_chrom$end_plus_CI_group/10^6, breaks=0:60)

  plot(NULL, xlim=c(0, max(consensus_df_chrom$end_plus_CI_group)), ylim=c(0, nrow(consensus_df_chrom)),
     main=chromosomes[i])
  arrows(x0=(consensus_df_chrom$start_plus_CI_group+
               consensus_df_chrom$start_minus_CI_group)/2, 
         x1=(consensus_df_chrom$end_plus_CI_group+
               consensus_df_chrom$end_minus_CI_group)/2,
         y0=(1:nrow(consensus_df_chrom)), 
         lwd=consensus_df_chrom$n_ind_group/20, col="red")
  
  
  hist(consensus_pass$start_minus_CI_group/10^6, breaks=0:60)
  hist(consensus_pass$end_plus_CI_group/10^6, breaks=0:60)
            
  plot(NULL, xlim=c(0, max(consensus_df_chrom$end_plus_CI_group)), ylim=c(0, npass),
     main=chromosomes[i])
  arrows(x0=(consensus_pass$start_plus_CI_group+
               consensus_pass$start_minus_CI_group)/2, 
       x1=(consensus_pass$end_plus_CI_group+
             consensus_pass$end_minus_CI_group)/2,
       y0=(1:npass), lwd=consensus_pass$n_ind_group/20+1, col="red")
dev.off()
 
}# end loop through chromosomes

sum(is.na(inv$ConsensusID)) #sanity check

# Write tables
write.csv(consensus_df,"../data_outputs/Consensus_List_unique.csv")
dim(consensus_df)
hist(consensus_df$ConsensusID, breaks=seq(0,30000, by=1000))

dim(inv)
dim(inv_new_df)
hist(inv_new_df$ConsensusID, breaks=seq(0,30000, by=1000))
write.csv(inv_new_df,"../data_outputs/INV-dataset-ConsensusID-unfiltered.csv")


### Filtered individuals DF ####
dim(inv_new_df)
ind_filtered <- merge(inv_new_df[,c("INDIVIDUAL","ConsensusID")], consensus_df, all.x=TRUE)
dim(ind_filtered) # same dimensions
head(ind_filtered)

ind_filtered <- ind_filtered[which(ind_filtered$pass),]
dim(ind_filtered)

ind_filtered <- ind_filtered[order(ind_filtered$INDIVIDUAL,
                                   ind_filtered$ConsensusID),]
head(ind_filtered,50)

ind_filtered$dup <- duplicated(ind_filtered)
sum(ind_filtered$dup)

### Final individual filtered ####
ind_filtered_final <- ind_filtered[!ind_filtered$dup,]
dim(ind_filtered_final)


ind_filtered_final$pop.names <- substr(ind_filtered_final$INDIVIDUAL,0,4)
ind_filtered_final$samp_full <- substr(ind_filtered_final$INDIVIDUAL,6,10)
tail(ind_filtered_final) 
head(ind_filtered_final) 

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

identical(ind_filtered_final3)

ind_filtered_final3 <- ind_filtered_final3[order(ind_filtered_final3$INDIVIDUAL,
                                   ind_filtered_final3$ConsensusID),]

write.csv(ind_filtered_final3,"../data_outputs/INV-dataset-ConsensusID-Ind-Filtered.csv")


#number of inversions per individual ####
inv_per_ind <- table(ind_filtered_final3$INDIVIDUAL)
hist(inv_per_ind)
summary(as.numeric(inv_per_ind))

str(ind_filtered_final3)

## Loop through chrom and make individual level plots ####
for (i in 1:length(chromosomes)){
  
  print(chromosomes[i])
  
  ind_filtered_final3_chrom <- ind_filtered_final3[ind_filtered_final3$chrom==chromosomes[i],]
  ind_filtered_final3_chrom <- ind_filtered_final3_chrom[order(ind_filtered_final3_chrom$chrom, 
                                                 ind_filtered_final3_chrom$start_minus_CI_group,
                                                 ind_filtered_final3_chrom$end_plus_CI_group,
                                                 ind_filtered_final3_chrom$pop),]
  head(ind_filtered_final3_chrom)
  
  pdf(paste0("../figures/",chromosomes[i],"ConsensusCalls_IndLevel.pdf"), height=20, width=10)
  
  plot(NULL, xlim=c(0, max(ind_filtered_final3_chrom$end_plus_CI_group)), 
       ylim=c(0, nrow(ind_filtered_final3_chrom)),
       main=chromosomes[i])
  arrows(x0=(ind_filtered_final3_chrom$start_plus_CI_group+
               ind_filtered_final3_chrom$start_minus_CI_group)/2, 
         x1=(ind_filtered_final3_chrom$end_plus_CI_group+
               ind_filtered_final3_chrom$end_minus_CI_group)/2,
         y0=(1:nrow(ind_filtered_final3_chrom)), 
         col=as.character(ind_filtered_final3_chrom$colorEcotypes))

  dev.off()
  
}# end loop through chromosomes
