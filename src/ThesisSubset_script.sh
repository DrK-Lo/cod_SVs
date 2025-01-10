#!/bin/bash
#SBATCH --job-name=Thesis_Subset
#SBATCH --mail-user=curtis.lei@northeastern.edu
#SBATCH --partition=lotterhos
#SBATCH --mem=35Gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-30%3
#SBATCH --output=/work/lotterhos/2020_CodGenomes_Inversions/2024Cod_Inversions/slurm_log/Thesis_Subset%j.out
#SBATCH --error=/work/lotterhos/2020_CodGenomes_Inversions/2024Cod_Inversions/slurm_log/Thesis_Subset%j.err

#Define variables
Cod_ID=`sed -n ${SLURM_ARRAY_TASK_ID}p /work/lotterhos/2020_CodGenomes_Inversions/2024Cod_Inversions/src/Cod_ID_ThesisSubset_Variables.txt`
REFERENCE="GCF_902167405.1_gadMor3.0_genomic"

## Load the necessary modules
module load singularity/3.10.3
module load samtools/1.19.2
echo "module load complete"

## Sort and index the BAM file
cd /home/curtis.lei
cp /work/lotterhos/2020_CodGenomes/labeled_bam_Out/${Cod_ID}.f.rg.bam /home/curtis.lei
samtools sort -m 40G ${Cod_ID}.f.rg.bam > ${Cod_ID}.sorted.f.rg.bam
samtools index ${Cod_ID}.sorted.f.rg.bam
echo "sort and index complete"

## Run DELLY on the sorted and indexed file
singularity run /home/curtis.lei/Containers/delly_v1.2.6.sif delly call -g ${REFERENCE}.fna ${Cod_ID}.sorted.f.rg.bam > /work/lotterhos/2020_CodGenomes_Inversions/outputs/delly/delly.${Cod_ID}.vcf
echo "delly complete"

## Run GRIDSS on generated  VCF
cp /work/lotterhos/2020_CodGenomes_Inversions/outputs/delly/delly.${Cod_ID}.vcf /home/curtis.lei
singularity run /home/curtis.lei/Containers/gridss_2.13.2.sif gridss_extract_overlapping_fragments --targetvcf delly.${Cod_ID}.vcf ${Cod_ID}.sorted.f.rg.bam
singularity run /home/curtis.lei/Containers/gridss_2.13.2.sif gridss -r ${REFERENCE}.fna -o gridssDelly.${Cod_ID}.vcf ${Cod_ID}.sorted.f.rg.bam.targeted.bam
echo "GRIDSS complete"

## Manage files in my home directory
mv gridssDelly.${Cod_ID}.vcf /work/lotterhos/2020_CodGenomes_Inversions/outputs/gridss_vcf
rm -r ${Cod_ID}.f.rg.bam ${Cod_ID}.sorted.f.rg.bam ${Cod_ID}.sorted.f.rg.bam.bai ${Cod_ID}.sorted.f.rg.bam.targeted.bam gridssDelly.${Cod_ID}.vcf.assembly.bam ${Cod_ID}.sorted.f.rg.bam.targeted.bam.gridss.working gridssDelly.${Cod_ID}.vcf.assembly.bam.gridss.working gridssDelly.${Cod_ID}.vcf.gridss.working delly.${Cod_ID}.vcf gridssDelly.${Cod_ID}.vcf.idx
echo "Files removed"
