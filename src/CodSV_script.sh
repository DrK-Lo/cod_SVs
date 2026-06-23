#!/bin/bash
#SBATCH --job-name=CodSV
#SBATCH --mail-user=curtis.lei@northeastern.edu
#SBATCH --partition=lotterhos
#SBATCH --mem=35Gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-296%5
#SBATCH --output=/projects/lotterhos/2020_CodGenomes_Inversions/2024Cod_Inversions/slurm_log/CodSV%j.out
#SBATCH --error=/projects/lotterhos/2020_CodGenomes_Inversions/2024Cod_Inversions/slurm_log/CodSV%j.err

#Load module and define variables
module load samtools/1.21
Cod_ID=`sed -n ${SLURM_ARRAY_TASK_ID}p /projects/lotterhos/2020_CodGenomes_Inversions/2024Cod_Inversions/src/Cod_ID_Variables.txt`
echo "module load and variables complete"

#Sort and index the BAM file
mkdir /scratch/curtis.lei/${Cod_ID} && cd /scratch/curtis.lei/${Cod_ID}
cp /projects/lotterhos/2020_CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna .
cp /projects/lotterhos/2020_CodGenomes/labeled_bam_Out/${Cod_ID}.f.rg.bam .
samtools sort -m 40G ${Cod_ID}.f.rg.bam > ${Cod_ID}.sorted.f.rg.bam
samtools index ${Cod_ID}.sorted.f.rg.bam
echo "sort and index complete"

#DELLY on .bam file
apptainer run /scratch/curtis.lei/delly_1.5.0.sif delly call -g GCF_902167405.1_gadMor3.0_genomic.fna ${Cod_ID}.sorted.f.rg.bam > delly.${Cod_ID}.vcf
echo "delly complete"

#GRIDSS on .bam file
apptainer run /scratch/curtis.lei/gridss_fedora.sif gridss -r GCF_902167405.1_gadMor3.0_genomic.fna -o gridss.${Cod_ID}.vcf ${Cod_ID}.sorted.f.rg.bam
eco "gridss complete"

#GRIDSS on delly.${Cod_ID}.vcf file
apptainer run /scratch/curtis.lei/gridss_fedora.sif gridss_extract_overlapping_fragments --targetvcf delly.${Cod_ID}.vcf ${Cod_ID}.sorted.f.rg.bam
apptainer run /scratch/curtis.lei/gridss_fedora.sif gridss -r GCF_902167405.1_gadMor3.0_genomic.fna -o gridssDelly.${Cod_ID}.vcf ${Cod_ID}.sorted.f.rg.bam.targeted.bam
echo "gridss_extract_overlapping_fragments complete"

## Manage files
mv delly.${Cod_ID}.vcf /projects/lotterhos/2020_CodGenomes_Inversions/outputs/delly
mv gridss.${Cod_ID}.vcf /projects/lotterhos/2020_CodGenomes_Inversions/outputs/gridss
mv gridssDelly.${Cod_ID}.vcf /projects/lotterhos/2020_CodGenomes_Inversions/outputs/gridss-delly
echo "Files moved"
rm -r /scratch/curtis.lei/${Cod_ID}/*
echo "Files removed"