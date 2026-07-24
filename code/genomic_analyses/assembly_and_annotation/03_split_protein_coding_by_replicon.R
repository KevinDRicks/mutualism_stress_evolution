#R script extracts all the predicted protein coding data from the prokka output and splits into seprate replicons
#splitting by replicon makes it more straight forward in downstream processes for mapping annotations to specific snps
#this is looped for each of the 28 ancestral strains

#load packages
library("Biostrings")

#read in prokka output
annotation_layout<-readLines("path/to/prokka/PROKKA_XXXXX.gff")
annotation_layout <- annotation_layout[-c(1:6)]

annot_layout_clean <- NULL
for(i in 1:length(annotation_layout)){
  working_string <- as.character(annotation_layout[i])
  cluster <- substr(working_string, start = 1, stop = 21)
  working_string <- sub("^[^,]*ID=", "", working_string)
  ID <- substr(working_string, start = 1, stop = 14)
  annot_layout_clean <- rbind(annot_layout_clean,c(cluster,ID))
}


protein_seq<-readAAStringSet("PROKKA_XXXXX.faa")

main_chromosome <- AAStringSet()
plasmidA <- AAStringSet()
plasmidB <- AAStringSet()
plasmidC <- AAStringSet()
plasmidD <- AAStringSet()

main_chromosome_ID <- NULL
plasmidA_ID <- NULL
plasmidB_ID <- NULL
plasmidC_ID <- NULL
plasmidD_ID <- NULL

#loops through all protein sequences and assigns them to a specific replicon
#number of plasmids are adjusted based on replicons
for(i in 1:length(protein_seq)){
  working_seq <- protein_seq[i]
  ID <- names(working_seq)
  ID_label <- substr(ID , start = 1, stop = 14)
  ID_protein <- substr(ID , start = 16, stop = nchar(ID))
  cluster_ID <- annot_layout_clean[annot_layout_clean[,2]==ID_label,1]
  
  if(cluster_ID=='cluster_001_consensus'){
    main_chromosome_ID <- c(main_chromosome_ID,ID_protein)
    main_chromosome <- c(main_chromosome,working_seq)
  }
  
  if(cluster_ID=='cluster_002_consensus'){
    plasmidA_ID <- c(plasmidA_ID,ID_protein)
    plasmidA <- c(plasmidA,working_seq)
  }
  
  if(cluster_ID=='cluster_003_consensus'){
    plasmidB_ID <- c(plasmidB_ID,ID_protein)
    plasmidB <- c(plasmidB,working_seq)
  }
  
  
  if(cluster_ID=='cluster_004_consensus'){
    plasmidC_ID <- c(plasmidC_ID,ID_protein)
    plasmidC <- c(plasmidC,working_seq)
  }
  
  if(cluster_ID=='cluster_005_consensus'){
    plasmidD_ID <- c(plasmidD_ID,ID_protein)
    plasmidD <- c(plasmidD,working_seq)
  } 
}


#outputs protein predictions
writeXStringSet(main_chromosome,'chromosome_AA.faa')
writeXStringSet(plasmidA,'plasmidA_AA.faa')
writeXStringSet(plasmidB,'plasmidB_AA.faa')
writeXStringSet(plasmidC,'plasmidC_AA.faa')
writeXStringSet(plasmidD,'plasmidD_AA.faa')

write.table(main_chromosome_ID,'chromosome_prokka.txt')
write.table(plasmidA_ID,'plasmidA_prokka.txt')
write.table(plasmidB_ID,'plasmidB_prokka.txt')
write.table(plasmidC_ID,'plasmidC_prokka.txt')
write.table(plasmidD_ID,'plasmidD_prokka.txt')