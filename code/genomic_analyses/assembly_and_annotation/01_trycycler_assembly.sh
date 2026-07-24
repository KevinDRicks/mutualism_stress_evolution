#code exectues the trycycler pipeline, building a full assembly for each Rhizobium strain
#this follow the standard pipeline for trycycler, which can be found in full detail here
#https://github.com/rrwick/Trycycler
#we use long read sequences (sequenced using PacBio Sequel II) in the assembly following followed by polishing
#raw sequencing files can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1294393
#note: Trycycler v0.5.5
#note: this script was initially run on Compute Canada server. Node time may be required for certain steps
#note: this script was looped for all our separate assemblies

#we created a virtual environment on the server termed "trycycler_assembly"
#this includes all libraries recommended by the pipeline (see trycycler github)

#activate trycycler virtual environemnt and set working environment
source ~/.virtualenvs/trycycler_assembly/bin/activate

#load libraries
module load StdEnv/2023
module load r/4.4.0
module load gcc/12.3
module load mash/2.3
module load muscle/3.8.1551

#copy in raw long read sequencing files
cat path/to/your/long_reads/*.fastq > raw_sequences.fastq
#pass long read sequences through standard filter, keeping only high quality reads
filtlong --min_length 1000 --keep_percent 95 raw_sequences.fastq > reads.fastq

#subsample sequences in order to build separate assemblies
trycycler subsample --reads reads.fastq --out_dir read_subsets 

#output directory for assemblies
mkdir assemblies

#with the subsample sequences, build separate assemblies using flye, raven, and minipolish
threads=16  #change as appropriate for your system

flye --nano-hq read_subsets/sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01
flye --nano-hq read_subsets/sample_02.fastq --threads "$threads" --out-dir assembly_02 && cp assembly_02/assembly.fasta assemblies/assembly_02.fasta && cp assembly_02/assembly_graph.gfa assemblies/assembly_02.gfa && rm -r assembly_02
flye --nano-hq read_subsets/sample_03.fastq --threads "$threads" --out-dir assembly_03 && cp assembly_03/assembly.fasta assemblies/assembly_03.fasta && cp assembly_03/assembly_graph.gfa assemblies/assembly_03.gfa && rm -r assembly_03
flye --nano-hq read_subsets/sample_04.fastq --threads "$threads" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_05.gfa read_subsets/sample_05.fastq > assemblies/assembly_05.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_06.gfa read_subsets/sample_06.fastq > assemblies/assembly_06.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_07.gfa read_subsets/sample_07.fastq > assemblies/assembly_07.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_08.gfa read_subsets/sample_08.fastq > assemblies/assembly_08.fasta

miniasm_and_minipolish.sh read_subsets/sample_09.fastq "$threads" > assemblies/assembly_09.gfa && any2fasta assemblies/assembly_09.gfa > assemblies/assembly_09.fasta
miniasm_and_minipolish.sh read_subsets/sample_10.fastq "$threads" > assemblies/assembly_10.gfa && any2fasta assemblies/assembly_10.gfa > assemblies/assembly_10.fasta
miniasm_and_minipolish.sh read_subsets/sample_11.fastq "$threads" > assemblies/assembly_11.gfa && any2fasta assemblies/assembly_11.gfa > assemblies/assembly_11.fasta
miniasm_and_minipolish.sh read_subsets/sample_12.fastq "$threads" > assemblies/assembly_12.gfa && any2fasta assemblies/assembly_12.gfa > assemblies/assembly_12.fasta


#cluster all the separate assemblies
trycycler cluster --assemblies assemblies/*.fasta --reads reads.fastq --out_dir trycycler
#manual examination of the clustering (as recommended by trycycler) was done to prune low quality assemblies
#5 clusters were identified, corresponding with  5 replicons
#note, this cluster number depended on the strain. Our ancestral strains had varying numbers of plasmids
#and consequently, each assembly had varying cluster content that is adjusted in each of these assemblies

#assemblies within each cluster were reconcilled
trycycler reconcile --reads reads.fastq --cluster_dir trycycler/cluster_001 --max_add_seq 20000
trycycler reconcile --reads reads.fastq --cluster_dir trycycler/cluster_002
trycycler reconcile --reads reads.fastq --cluster_dir trycycler/cluster_003
trycycler reconcile --reads reads.fastq --cluster_dir trycycler/cluster_004
trycycler reconcile --reads reads.fastq --cluster_dir trycycler/cluster_005

#sequence alignments are run against the reconciled clusters
trycycler msa --cluster_dir trycycler/cluster_001
trycycler msa --cluster_dir trycycler/cluster_002
trycycler msa --cluster_dir trycycler/cluster_003
trycycler msa --cluster_dir trycycler/cluster_004
trycycler msa --cluster_dir trycycler/cluster_005


#partition the long read sequences between these clusters 
trycycler partition --reads reads.fastq --cluster_dirs trycycler/cluster_*

#Generate a consensus assembly for each cluster
trycycler consensus --cluster_dir trycycler/cluster_001
trycycler consensus --cluster_dir trycycler/cluster_002
trycycler consensus --cluster_dir trycycler/cluster_003
trycycler consensus --cluster_dir trycycler/cluster_004
trycycler consensus --cluster_dir trycycler/cluster_005

#polish the consensus assembly with long reads using medaka
#note, may require node time
module load parasail
module load bcftools/1.19
module load samtools/1.20
module load htslib/1.19
module load minimap2

medaka_consensus -i trycycler/cluster_001/4_reads.fastq -d trycycler/cluster_001/7_final_consensus.fasta -o trycycler/cluster_001/medaka -m r941_min_sup_g507 -t 12
medaka_consensus -i trycycler/cluster_002/4_reads.fastq -d trycycler/cluster_002/7_final_consensus.fasta -o trycycler/cluster_002/medaka -m r941_min_sup_g507 -t 12
medaka_consensus -i trycycler/cluster_003/4_reads.fastq -d trycycler/cluster_003/7_final_consensus.fasta -o trycycler/cluster_003/medaka -m r941_min_sup_g507 -t 12
medaka_consensus -i trycycler/cluster_004/4_reads.fastq -d trycycler/cluster_004/7_final_consensus.fasta -o trycycler/cluster_004/medaka -m r941_min_sup_g507 -t 12
medaka_consensus -i trycycler/cluster_005/4_reads.fastq -d trycycler/cluster_005/7_final_consensus.fasta -o trycycler/cluster_005/medaka -m r941_min_sup_g507 -t 12

