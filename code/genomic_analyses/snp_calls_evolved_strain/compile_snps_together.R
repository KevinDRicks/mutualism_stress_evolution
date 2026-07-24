#code takes the outputs from the spine, gatk and breseq pipeline and
#transforms them into a single harmonized snp table, associating each snp with the proper treatments
#we additionally generate summary statistic here 
#output file is the snp_summary.csv, which is used in downstream analysis


library(readr)

setwd('G:/My Drive/Work/UIUC/Projects/Rhizobia/Drought_evolution/Data/2025_update/genomics/snps')
list.files()
spine_snp <- read.csv('spine_out.csv')
breseq_snp <- read.table('breseq_out.txt',header=T)
gatk_snp <- read.table('gatk_out.txt',header=T)
mapper <- read.csv('evolved_strain_identity.csv')
ancestor.mapper <- read.csv('ancestor_label.csv')


#massaging gatk breseq and spine files into the same format 
gatk_snp$File <- gsub(".new.vcf","",gatk_snp$File)
colnames(gatk_snp) <- c('Evolved','Cluster','Position','Reference','Variant')
gatk_snp$Sample <- gsub("KR24-","",gatk_snp$Evolved)
gatk_snp$Sample <- sub("\\..*", "", gatk_snp$Sample)
gatk_snp$Sample <- paste('Sample',gatk_snp$Sample,sep='')
gatk_snp$Ancestor <- NA
for(i in 1:nrow(gatk_snp)){
  gatk_snp$Ancestor[i] <- mapper[mapper$Sequencing_name==gatk_snp$Sample[i],]$Ancestor
}
colnames(breseq_snp)
breseq_snp <- data.frame(Evolved=breseq_snp$Evolved,Cluster=breseq_snp$Cluster,Position=breseq_snp$Position,
                         Reference=breseq_snp$Ancestral_state,Variant=breseq_snp$Evolved_state,Sample=NA,
                         Ancestor= breseq_snp$Ancestor)
spine_snp <- data.frame(Evolved=spine_snp$Evolved,Cluster=spine_snp$Cluster,Position=spine_snp$Position,
                        Reference=spine_snp$Ancestral_state,Variant=spine_snp$Evolved_state,Sample=NA,
                        Ancestor= spine_snp$Ancestor)


spine_snp$Sample <- paste('Sample',sub("\\..*", "", gsub("KR24-","",spine_snp$Evolved)),sep='')
breseq_snp$Sample <- paste('Sample',sub("\\..*", "", gsub("KR24-","",breseq_snp$Evolved)),sep='')



spine_snp$Ancestor <- gsub("Rht_","",spine_snp$Ancestor)
breseq_snp$Ancestor <- gsub("Rht_","",breseq_snp$Ancestor)
spine_snp$Ancestor <- gsub("_C","",spine_snp$Ancestor)
breseq_snp$Ancestor <- gsub("_C","",breseq_snp$Ancestor)
spine_snp$Ancestor <- gsub("_N","",spine_snp$Ancestor)
breseq_snp$Ancestor <- gsub("_N","",breseq_snp$Ancestor)

length(unique(spine_snp$Evolved))
length(unique(breseq_snp$Evolved))
length(unique(gatk_snp$Evolved))
length(unique(c(spine_snp$Evolved,breseq_snp$Evolved,gatk_snp$Evolved)))

breseq_snp_short <- breseq_snp[nchar(breseq_snp$Reference)<50,]
breseq_snp_long <- breseq_snp[nchar(breseq_snp$Reference)>50,]
length(breseq_snp_long)

breseq_snp_long[,-4]
breseq_snp_long[breseq_snp_long$Ancestor=='262262_plasmid1,length=1227174>',]$Ancestor <- '262'
breseq_snp_short <- rbind(breseq_snp_short,breseq_snp_long)

nchar(breseq_snp_long[,4])


merged_snp <- spine_snp
dim(merged_snp)
for(i in 1:nrow(breseq_snp_short)){
  working.smp <- breseq_snp_short[i,]
  if(sum(working.smp$Evolved==merged_snp$Evolved & abs(working.smp$Position-merged_snp$Position)<5)==0){
    merged_snp <- rbind(merged_snp,working.smp)
  }
}
dim(merged_snp)#416 snps



for(i in 1:nrow(gatk_snp)){
  working.smp <- gatk_snp[i,]
  if(sum(working.smp$Evolved==merged_snp$Evolved & abs(working.smp$Position-merged_snp$Position)<5)==0){
    merged_snp <- rbind(merged_snp,working.smp)
  }
}
dim(merged_snp)#416 snps



#gatk frequently called single deletions where breseq called whole plasmid deletions
#clean those up and remove from sample 
for(i in 1:nrow(breseq_snp_long)){
  pos <- which(merged_snp$Evolved==breseq_snp_long[i,]$Evolved & merged_snp$Cluster==breseq_snp_long[i,]$Cluster & 
                 merged_snp$Position>breseq_snp_long[i,]$Position & merged_snp$Position<breseq_snp_long[i,]$EndPosition)
  if(length(pos)>0){
    merged_snp <- merged_snp[-pos,]
  }
}
dim(merged_snp)#288 snps

sort(table(merged_snp$Evolved))



#make sure that there are no duplicates
new.subset <- merged_snp[1,]
for(i in 2:nrow(merged_snp)){
  if(sum(new.subset$Strain==merged_snp$Strain[i] & new.subset$Position==merged_snp$Position[i])==0){
    new.subset <- rbind(new.subset,merged_snp[i,])
  }
}
merged_snp <- new.subset 
dim(merged_snp)


merged_snp_exclude <- merged_snp#[!(merged_snp$Evolved=='KR24-20.19_19' | merged_snp$Evolved=='KR24-34.33_33'), ]
dim(merged_snp_exclude)
table(merged_snp_exclude$Cluster)

merged_snp_exclude$Strain <- NA
for(i in 1:nrow(merged_snp_exclude)){
  merged_snp_exclude$Strain[i] <- mapper[mapper$Sequencing_name==merged_snp_exclude$Sample[i],]$Greenhouse_name
}


write.csv(merged_snp_exclude,"G:/My Drive/Work/UIUC/Projects/Rhizobia/Drought_evolution/Data/2025_update/genomics/snp.csv",row.names = F)
