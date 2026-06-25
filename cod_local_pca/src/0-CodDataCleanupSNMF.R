library(bigstatsr)
library(bigsnpr)
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

### Impute genotypes
write.geno(G, "~/Desktop/codGenotypeData/merged.f.99ind.MAF05_imputed.geno")
project.snmf = snmf("~/Desktop/codGenotypeData/merged.f.99ind.MAF05_imputed.geno", K = 1:5, 
                    repetitions=5,
                    entropy = TRUE, project = "new")

project.snmf = load.snmfProject("~/Desktop/codGenotypeData/merged.f.99ind.MAF05_imputed.snmfProject")

plot(project.snmf, lwd = 5, col = "red", pch=1)
ce = cross.entropy(project.snmf, K = 3)
best_run = which.min(ce)
best_run
# 4. Impute missing values
impute(project.snmf, 
                    "~/Desktop/codGenotypeData/merged.f.99ind.MAF05_imputed.geno", 
       method = "genotype", 
       K = 3, 
       run = best_run)
