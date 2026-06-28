
### Starting code
SV_start = 19
SV_end = 30
PCA_start = 9
PCA_end = 22
(SV_size  =  SV_end - SV_start)
(PCA_size = PCA_end - PCA_start)

# Does at least 70% of SV range overlap with PCA range
SV_list = SV_start:SV_end
PCA_list = PCA_start:PCA_end
Overlap = intersect(SV_list, PCA_list) # get the intersection of the two ranges
(PercentOverlapSV = length(Overlap)/length(SV_list)) # get the percent of the SV range that is in the overlap
(DoesOverlap = PercentOverlapSV >= 0.7) #logical used to evaluate final overlap
# The 0.7 is a threshold that at least 70% of the SV call has to overlap with the PCA range
# (change threshold as needed)

# Is SV size within 0.5 to 1.5 the size of PCA size
(SizeRatio = SV_size/PCA_size)
(DoesSizeMatch = SizeRatio >= 0.5 & SizeRatio <= 1.5)

# Final assessment
(IsItAMatch <- DoesOverlap & DoesSizeMatch)


evaluateAgreementSV_PCA <- function{PCA_genotype, anySV_present, SV_start, SV_end,
  PCA_start, PCA_end}

# Example values
PCA_genotype = 0
SV_present = TRUE

SV_start = -5
SV_end = 50
PCA_start = 0
PCA_end = 100

if (PCA_genotype == 0) and (SV_present == FALSE){
  IsItAMatch = "PCA absent, SV absent"
  break
}

if (PCA_genotype == 0) and (SV_present == TRUE){
  IsItAMatch = "PCA absent, SV present"
  break
}

if (PCA_genotype > 0) and (SV_present == FALSE){
  IsItAMatch = "PCA present, SV absent"
  break
}

if (PCA_genotype > 0) and (SV_present == TRUE){
  
  (SV_size  =  SV_end - SV_start)
  (PCA_size = PCA_end - PCA_start)
  
  # Does at least 70% of SV range overlap with PCA range
  SV_list = SV_start:SV_end
  PCA_list = PCA_start:PCA_end
  Overlap = intersect(SV_list, PCA_list) # get the intersection of the two ranges
  (PercentOverlapSV = length(Overlap)/length(SV_list)) # get the percent of the SV range that is in the overlap
  (DoesOverlap = PercentOverlapSV >= 0.7) #logical used to evaluate final overlap 
   # The 0.7 is a threshold that at least 70% of the SV call has to overlap with the PCA range 
    # (change threshold as needed)
  
  # Is SV size within 0.5 to 1.5 the size of PCA size
  (SizeRatio = SV_size/PCA_size)
  (DoesSizeMatch = SizeRatio >= 0.5 & SizeRatio <= 1.5)
  
  
  DoesOverlap & DoesSizeMatch
  # Final assessment
   if (DoesOverlap & DoesSizeMatch){
     IsItAMatch =  "PCA present, SV present"
   }else{
     IsItAMatch = "PCA present, SV absent"
   }
 
  return(IsItAMatch) 
}

#####

#Lei Code
Size_Pass <- SV_start <= PCA_size+(PCA_size*0.5) & SV_start <= PCA_size+(PCA_size*0.5)
(Location_Pass <- abs((SV_start-SV_end)/(pmin(SV_start,PCA_start)-pmax(SV_end,PCA_end)))>=0.7)
ID_PASS <- Size_Pass=="TRUE" & Location_Pass=="TRUE"