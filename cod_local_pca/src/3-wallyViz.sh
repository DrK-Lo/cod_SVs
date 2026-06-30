

# get off the login node
srun -p short -N 1 --pty /bin/bash


# pull the docker container (DO ONCE, NO NEED TO RERUN)
#cd /projects/lotterhos/Containers
#ls
# export APPTAINER_CACHEDIR=$(pwd)/cache 
# export APPTAINER_TMPDIR=$(pwd)/tmp
#apptainer pull wally.sif docker://trausch/wally
#ls

# RUN THE APPTAINER
cd /projects/lotterhos/Containers
apptainer run -B $(cd ../ && pwd):/projects/lotterhos wally.sif

# run wally with split
path_out="/projects/lotterhos/2020_CodGenomes/labeled_bam_Out/ReadVisualizations/"
cd ${path_out}

# example of working wally code (uncomment to run)
#wally region -p -r NC_044048.1:11250000-11500000:myplot -g ../../Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna ../indexedSorted/Pop1_16216.f.rg.bam_sorted.bam
#wally region -p -r NC_044048.1:11250000-11500000:/projects/lotterhos/2020_CodGenomes/labeled_bam_Out/ReadVisualizations/NC_044048.1_11250000-11500000/NC_044048.1_11250000-11500000_Pop1_16216.f.rg.bam -g ../../Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna ../indexedSorted/Pop1_16216.f.rg.bam_sorted.bam 
#wally region -p -r NC_044048.1:11250000-11500000 -g ../../Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna ../indexedSorted/Pop1_16216.f.rg.bam_sorted.bam 

genome=../../Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna
echo $genome

# input the region and create a directory to host the files
#split=1
#myregion=${chr}:11250000-11600000 #example with one window

### This is the section to edit for where you want to see the mapped reads
chr=NC_044048.1
split=4
myregion=${chr}:11250000-11500000,${chr}:13200000-13500000,${chr}:20500000-21000000,${chr}:28000000-29000000 #xample with split window
myplot1=$(echo $myregion | tr ':,' '_')
echo $myplot1
mkdir ${myplot1}

# loop through all the bams and create a visualizations
for i in {1..296}; do

	ID=$(sed -n "${i}p" ../bamlist.txt)
	#echo $ID

	mybam=../indexedSorted/${ID}_sorted.bam
	#echo ${mybam}

	myplot2="${path_out}${myplot1}/${myplot1}_${ID}"
	#echo $myplot2

	wally region -s ${split} -p -r ${myregion}:${myplot2} -g ${genome} ${mybam} -x 3000
done



###################################
###################################
# Index and sort a bam ####
## SortINdexBams.sh ####

#!/bin/bash

#SBATCH --job-name=CodBamIndexSort
#SBATCH --mem=20Gb
#SBATCH --mail-user=k.lotterhos@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-296%50
#SBATCH --output=/projects/lotterhos/2020_CodGenomes/labeled_bam_Out/indexedSorted/CodBamIndexSort.out
#SBATCH --error=/projects/lotterhos/2020_CodGenomes/labeled_bam_Out/indexedSorted/CodBamIndexSort.err

module load samtools
cd /projects/lotterhos/2020_CodGenomes/labeled_bam_Out/indexedSorted

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../bamlist.txt)
echo $ID

samtools sort "../${ID}" -o ${ID}_sorted.bam
echo "samtools sort done"
samtools index ${ID}_sorted.bam
echo "samtools index done"

###################################
###################################

### [lotterhos@explorer-02 labeled_bam_Out]$ 

srun -p short -N 1 --pty /bin/bash
cd /projects/lotterhos/2020_CodGenomes/labeled_bam_Out/indexedSorted
vim ../SortIndexBams.sh

sbatch ../SortIndexBams.sh


rm -f *.bam.tmp.*.bam
sbatch ../SortIndexBams.sh
squeue -u lotterhos




###################################
###################################