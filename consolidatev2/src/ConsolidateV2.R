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

# study CI ####
hist(inv$POS_PLUS_CI-inv$POS_MINUS_CI)
hist(inv$END_PLUS_CI-inv$END_MINUS_CI)
plot(inv$POS_PLUS_CI-inv$POS_MINUS_CI, inv$END_PLUS_CI-inv$END_MINUS_CI)
plot(inv$maxSize, inv$POS_PLUS_CI-inv$POS_MINUS_CI)



## Sort inv by chromsome,POS_PLUS_CI## Sort inv by chromsome, pos, end
myorder <- order(inv$CHROM, inv$POS_MINUS_CI, inv$END_PLUS_CI)
inv <- inv[myorder,]
head(inv)
inv$ConsensusID <- NA


## Inspect FILTER
levels(as.factor(inv$FILTER))
tapply(inv$FILTER, inv$FILTER, length)

### Being looping here, reset vars
consensus_df = NULL
ConsensusID_counter <- 1
overlap <- 0
inv$ConsensusID = NA

### Loop through chromosomes
for (i in 1:length(chromosomes)){
  
  inv_indices <- inv$CHROM==chromosomes[i]
  inv_chrom <- inv[inv_indices,]
  print(dim(inv_chrom))
  
  pdf(paste0("../figures/",chromosomes[i],"AllCalls.pdf"), height=20, width=10)
  plot(NULL, xlim=c(0, max(inv_chrom$END_PLUS_CI)), ylim=c(0, nrow(inv_chrom)))
  arrows(x0=(inv_chrom$POS_MINUS_CI), 
         x1=inv_chrom$END_PLUS_CI,
       y0=(1:nrow(inv_chrom)), col=(inv_chrom$PROGRAM_col), xlim=c(0,30*10^6),
       main=chromosomes[i], lwd=0.1)
  plot(NULL, xlim=c(0, max(inv_chrom$END_PLUS_CI)), ylim=c(0, nrow(inv_chrom)))
  arrows(x0=(inv_chrom$POS_MINUS_CI), 
         x1=inv_chrom$END_PLUS_CI,
           y0=(1:nrow(inv_chrom)), col=(inv_chrom$FILTER=="PASS")+1, xlim=c(0,30*10^6),
         main=chromosomes[i], lwd=0.1)
  dev.off()
  
  while(sum(is.na(inv_chrom$ConsensusID)) > 0){
    
    (index = min(which(is.na(inv_chrom$ConsensusID))))
    
    print(c(index,"of", nrow(inv_chrom)))
    
    (focal_start_minus_CI <- inv_chrom$POS_MINUS_CI[index])
    (focal_start_plus_CI <- inv_chrom$POS_PLUS_CI[index])
    (focal_end_plus_CI <- inv_chrom$END_PLUS_CI[index])
    (focal_end_minus_CI <- inv_chrom$END_MINUS_CI[index])
    (focal_size <- inv_chrom$maxSize[index])
    
    ## The first group is based on a focal chromosome
    group_1 <- which((inv_chrom$POS_MINUS_CI+inv_chrom$POS_PLUS_CI)/2<=focal_start_plus_CI+overlap*focal_size &
                    (inv_chrom$POS_MINUS_CI+inv_chrom$POS_PLUS_CI)/2>=focal_start_minus_CI-overlap*focal_size &
                    (inv_chrom$END_MINUS_CI+inv_chrom$END_PLUS_CI)/2<=focal_end_plus_CI+overlap*focal_size &
                    (inv_chrom$END_MINUS_CI+inv_chrom$END_PLUS_CI)/2>=focal_end_minus_CI-overlap*focal_size 
    )
    
    #inv_chrom[group_1,]
    dim(inv_chrom[group_1,])
    
    # The second group is based on the median breakpoint location of the first group
    (focal_start_minus_CI2 <- min(inv_chrom$POS_MINUS_CI[group_1], na.rm=TRUE))
    (focal_start_plus_CI2 <- max(inv_chrom$POS_PLUS_CI[group_1], na.rm=TRUE))
    (focal_end_plus_CI2 <- max(inv_chrom$END_PLUS_CI[group_1], na.rm=TRUE))
    (focal_end_minus_CI2 <- min(inv_chrom$END_MINUS_CI[group_1], na.rm=TRUE))
    (focal_size2 <- focal_end_plus_CI2 -focal_start_minus_CI2)
    group_2 <- which((inv_chrom$POS_MINUS_CI+inv_chrom$POS_PLUS_CI)/2<=focal_start_plus_CI2+overlap*focal_size2 &
                       (inv_chrom$POS_MINUS_CI+inv_chrom$POS_PLUS_CI)/2>=focal_start_minus_CI2-overlap*focal_size2&
                       (inv_chrom$END_MINUS_CI+inv_chrom$END_PLUS_CI)/2<=focal_end_plus_CI2+overlap*focal_size2 &
                       (inv_chrom$END_MINUS_CI+inv_chrom$END_PLUS_CI)/2>=focal_end_minus_CI2-overlap*focal_size2 
    )
    
    # Give them the same ID
    dim(inv_chrom[group_2,])
    inv_chrom$ConsensusID[group_2] <- ConsensusID_counter
    
    # Calc summary stats
    (n_calls_group <- length(group_2))
    (n_programs_group <- length(levels(as.factor(inv_chrom$PROGRAM[group_2]))))
    (n_PASS_group <- sum(inv_chrom$FILTER[group_2]=="PASS"))
    (n_ind_group <- length(levels(as.factor(inv_chrom$INDIVIDUAL[group_2]))))
    
    (start_plus_CI_group = max(inv_chrom$POS_PLUS_CI[group_2], na.rm=TRUE))
    (start_minus_CI_group = min(inv_chrom$POS_MINUS_CI[group_2], na.rm=TRUE))
    (end_plus_CI_group = max(inv_chrom$END_PLUS_CI[group_2], na.rm=TRUE))
    (end_minus_CI_group = min(inv_chrom$END_MINUS_CI[group_2], na.rm=TRUE))
    
    # Save summary stats
    consensus_df_group <- data.frame(chrom = chromosomes[1],
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
    
    inv$ConsensusID[inv_indices] <- ConsensusID_counter
    ConsensusID_counter <- ConsensusID_counter + 1
  } #end while loop
    
  
  
}

head(consensus_df)
dim(consensus_df)

summary(consensus_df$start_minus_CI_group/10^6)
summary(consensus_df$start_plus_CI_group/10^6)
summary(consensus_df$start_plus_CI_group/10^6)


consensus_df$pass <- consensus_df$n_programs_group > 0 &
                      consensus_df$n_PASS_group >= 0 &
                      consensus_df$n_ind_group > 15

tapply(consensus_df$chrom, consensus_df$pass, length)

mypass <- which(consensus_df$pass)
npass <- length(mypass)


consensus_pass <- consensus_df[mypass,]
dim(consensus_pass)
consensus_pass <- consensus_df[order(consensus_pass$chrom, 
                                     consensus_pass$start_minus_CI_group,
                                     consensus_pass$end_plus_CI_group),]

pdf(paste0("../figures/",chromosomes[i],"ConsensusCalls.pdf"), height=20, width=10)
  plot(NULL, xlim=c(0, max(inv_chrom$END_PLUS_CI)), ylim=c(0, npass),
     main=chromosomes[i])
  arrows(x0=(consensus_pass$start_plus_CI_group+consensus_pass$start_minus_CI_group)/2, 
       x1=(consensus_pass$end_plus_CI_group+consensus_pass$end_minus_CI_group)/2,
       y0=(1:npass), lwd=consensus_pass$n_ind_group/20, col="red")
dev.off()


inv2 <- merge(inv, )





