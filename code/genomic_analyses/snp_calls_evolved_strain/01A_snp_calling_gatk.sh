#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --output=output_gatk.txt  
#SBATCH --error=error_gatk.txt
#SBATCH --mem=96G

#this script is one of multiple methods used to detect introduced mutations in our evolved strains
#in this script, raw sequences from evolved strains are passed through the gatk pipeline 
#and compared against the high quality ancestor assembly that have previously been built
#raw long reads for can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1294393/
#ancestral assemblies can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA310138/

#this script must by run for each evolved strain (modified to run in a loop)

########################################################################################################v
#Reference Genome Preparation 
#note: You only need to run this once per reference genome. If the ancestral assembly
is already indexed, you can comment this section out
########################################################################################################v
module load StdEnv/2023
module load bwa/0.7.18
ancestor_path="/path/to/ancestor/assembly/anc_assembly.fasta"
ancestor_dict="output/path/for/genome/dictionary"
bwa index $ancestor_path

#to run properly,  this required to separate StdEnv in the computecanada space 
#2023 above, and then 2020 below
#other clusters will likely have other requirements
module load StdEnv/2020
module load samtools/1.17 
module load picard/2.26.3
samtools faidx $ancestor_path
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary REFERENCE=$ancestor_path OUTPUT=$ancestor_dict


########################################################################################################v
#merge cleaned, and quality controlled short reads
########################################################################################################v
header="/home/kricks/scratch/saline_evolution_update/cleaned_reads/"
target_evolve1="/path/to/long_reads/[STRAIN_ID]_R1_longread_001.fastq.gz"
target_evolve2="/path/to/long_reads/[STRAIN_ID]_R2_longread_001.fastq.gz"
out_merged_gz="[STRAIN_ID]_merge.fastq.gz"
out_merged_unzipped="[STRAIN_ID]_merge.fastq"

cat $target_evolve1 $target_evolve2 > $out_merged_gz
gunzip $out_merged_gz


########################################################################################################v
#align reads to ancestor reference & convert to bam
########################################################################################################v
module load StdEnv/2023
module load bwa/0.7.18
module load samtools/1.20

bwa mem -t 4 -M -R "@RG\tID:[STRAIN_ID]\tSM:[STRAIN_ID]" $ancestor_path $out_merged_unzipped | samtools view -huS -o "[STRAIN_ID].bam" 


########################################################################################################v
#reorder using picard
########################################################################################################v
module load picard/3.1.0
module load java/21.0.1 

java -jar $EBROOTPICARD/picard.jar ReorderSam \
  R=$ancestor_path \
  SEQUENCE_DICTIONARY=$ancestor_dict \
  I="[STRAIN_ID].bam" \
  O="[STRAIN_ID].reorder.bam"


########################################################################################################v
#add read groups
########################################################################################################v
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
  -I "[STRAIN_ID].reorder.bam" \
  -O "[STRAIN_ID].new_rg.bam" \
  -ID "[STRAIN_ID]" \
  -LB saline_evo \
  -PL Illumina \
  -PU 1 \
  -SM "[STRAIN_ID]"


########################################################################################################v
#sort the input BAM file by coordinate
########################################################################################################v
java -jar $EBROOTPICARD/picard.jar SortSam \
  I="[STRAIN_ID].new_rg.bam" \
  O="[STRAIN_ID].coordinate_sorted.bam" \
  SO=coordinate

########################################################################################################v
#identify duplicate reads
########################################################################################################v
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  -I "[STRAIN_ID].coordinate_sorted.bam" \
  -O "[STRAIN_ID].marked_duplicates.bam" \
  -M "[STRAIN_ID].marked_dup_metrics.txt"

########################################################################################################v
#index the newly marked BAM file
########################################################################################################v
java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
  -I "[STRAIN_ID].marked_duplicates.bam"

########################################################################################################v
#generate stats
########################################################################################################v
samtools stats "[STRAIN_ID].marked_duplicates.bam" > "[STRAIN_ID].stats"


########################################################################################################v
#GATK haplotypeCaller
########################################################################################################v
GATK_JAR="/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar"

java -jar $GATK_JAR HaplotypeCaller \
  -R $ancestor_path \
  -I "[STRAIN_ID].marked_duplicates.bam" \
  --dont-use-soft-clipped-bases TRUE \
  -ploidy 1 \
  -O "[STRAIN_ID].g.vcf.gz" \
  -ERC GVCF

########################################################################################################v
#combine
########################################################################################################v
#once run for all strains, can combine files together
ls *.g.vcf.gz > gVCFs.list
java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar GenotypeGVCFs -R $ancestor_path -V ref.vcf -ploidy 1 -O new.vcf -stand-call-conf 30 > genotype_gVCFs.out 
java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar VariantsToTable -V new.vcf -R $ancestor_path -F CHROM -F POS -F REF -F ALT -F QUAL -F AF -F ANN -F DP -GF GT -O gatk_snps.table