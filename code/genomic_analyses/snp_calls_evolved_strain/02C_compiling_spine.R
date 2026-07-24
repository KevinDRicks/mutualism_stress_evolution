#following looping the spine script, an output file will be generated for every derived strain
#we loop through those an compile them into one large summary of the snps observed across the experiment


setwd('/path/to/spine/output')

files <- list.files()


evolved_df <- read.table('/path/to/evolved/strain/list.txt')
ancestral_df <- read.table('path/to/ancestral/strain/list.txt')


summary_df <- NULL
detail_df <- NULL
for(i in 1:length(files)){
  
  
  
  tmp_snp<-readLines(paste(files[i],'/NUCMER/b_ancestor.fasta_core.snps',sep=''))
  
  evolve_strain_tmp <- files[i]
  anc <- ancestral_df[evolved_df==evolve_strain_tmp,]
  
  
  working_df <- NULL
  if(length(tmp_snp)==5){
    vec_out <- c(anc,evolve_strain_tmp,0,0,0)
    names(vec_out) <- c('Ancestor_strain','Evolved','Num_unique_SNP','Num_unique_IN','Num_unique_DEL')  
    summary_df <- rbind(summary_df,vec_out)
  }
  
  if(length(tmp_snp)>5){

    working_df <- NULL
    for(j in 6:length(tmp_snp)){
      working_file <- tmp_snp[j]
      working_file <- unlist(strsplit(working_file, split = " "))
      working_file <- working_file[working_file!=""]
      
      output.d <- c(working_file[3],working_file[2],working_file[4])
      output.d[4] <- ifelse(output.d[1]!="." & output.d[2]!="." ,'y','n')
      output.d[5] <- ifelse(output.d[1]=='.','y','n')
      output.d[6] <- ifelse(output.d[2]=='.','y','n')
      cluster <- working_file[length(working_file)]
      cluster <- unlist(strsplit(cluster,split="\t"))[2]
      output.d <- c(anc,evolve_strain_tmp,cluster,output.d)
      working_df <- rbind(working_df,output.d)
    }
    
    
    
    vec_out <- c(anc,files[i],NA,NA,NA)
    names(vec_out) <- c('Ancestor_strain','Evolved','Num_unique_SNP','Num_unique_IN','Num_unique_DEL')  
    vec_out[3] <- sum(working_df[,7]=='y')
    vec_out[4] <- sum(working_df[,8]=='y')
    vec_out[5] <- sum(working_df[,9]=='y')
    summary_df <- rbind(summary_df,vec_out)
    
    
    
    colnames(working_df) <- c('Ancestor_strain','Evolved','Cluster','Ancestral_state','Evolved_state','Position','SNP','Insertion','Deletion')
    rownames(working_df) <- NULL
    working_df <- data.frame(working_df)
    working_df$Position <- as.numeric(working_df$Position)
    detail_df <- rbind(detail_df,working_df)
  }
}


#renaming to be correct cluster name
rename.df <- data.frame(detail_df,cluster_len=NA)
for(i in 1:nrow(rename.df)){
  setwd(paste('spine_output/',rename.df$Evolved[i],"/SPINE/",sep=""))
  working.d <- read.table('output.b_ancestor.fasta.core_coords.txt',header = T)
  rename.df$cluster_len[i] <-  working.d[working.d$out_seq_id==rename.df$Cluster[i],]$contig_length
  rename.df$Cluster[i] <-  working.d[working.d$out_seq_id==rename.df$Cluster[i],]$contig_id  
}
write.csv(rename.df,'spine_out.csv',row.names=F)
write.csv(summary_df,'spine_summary.csv',row.names=F)
