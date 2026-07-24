#following looping the breseq script, an output file will be generated for every derived strain
#we loop through those an compile them into one large summary of the snps observed across the experiment

setwd('/path/to/breseq/output')
files <- list.files()

files <- files[files!="breseq_out.txt"]
files <- files[files!="breseq_summary.txt"]
summary_df <- NULL
detail_df <- NULL
for(i in 1:length(files)){
  target_directory <- paste(files[i],"/output/output.vcf",sep="")
  
  if(file.exists(file = target_directory)){
    print(i)
    tmp_vcf<-readLines(target_directory)
    
    end_pos <- which(tmp_vcf=="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    working_df <- NULL
    if(end_pos==length(tmp_vcf)){
      anc <- tmp_vcf[5]
      anc<-gsub('##contig=<ID=','',anc)
      anc<-gsub('cluster_.*','',anc)
      vec_out <- c(anc,files[i],0,0,0)
      names(vec_out) <- c('Ancestor_strain','Evolved','Num_unique_SNP','Num_unique_IN','Num_unique_DEL')  
      summary_df <- rbind(summary_df,vec_out)
    }
    
    if(end_pos<length(tmp_vcf)){
      anc <- tmp_vcf[5]
      anc<-gsub('##contig=<ID=','',anc)
      anc<-gsub('cluster_.*','',anc)
      
      working_df <- NULL
      for(j in (end_pos+1):length(tmp_vcf)){
        working_file <- tmp_vcf[j]
        working_file <- unlist(strsplit(working_file, split = "\t"))
        
        output.d <- c(working_file[1],working_file[2],working_file[4],working_file[5])
        output.d[5] <- ifelse(nchar(output.d[3])==nchar(output.d[4]) & 
                              output.d[3]!=output.d[4],'y','n')
        output.d[6] <- ifelse(nchar(output.d[3])>nchar(output.d[4]),'y','n')
        output.d[7] <- ifelse(nchar(output.d[3])<nchar(output.d[4]),'y','n')
        working_df <- rbind(working_df,output.d)
      }
      
      
      
      vec_out <- c(anc,files[i],NA,NA,NA)
      names(vec_out) <- c('Ancestor_strain','Evolved','Num_unique_SNP','Num_unique_IN','Num_unique_DEL')  
      vec_out[3] <- sum(working_df[,5]=='y')
      vec_out[4] <- sum(working_df[,7]=='y')
      vec_out[5] <- sum(working_df[,6]=='y')
      summary_df <- rbind(summary_df,vec_out)
      
      
      
      colnames(working_df) <- c('Cluster','Position','Ancestral_state','Evolved_state','SNP','Insertion','Deletion')
      rownames(working_df) <- NULL
      working_df <- data.frame(working_df)
      working_df$Position <- as.numeric(working_df$Position)
      working_df$Position_upstream <- working_df$Position+40
      working_df$Position_downstream <- working_df$Position-40
      working_df$Position_adj <- working_df$Position + nchar(working_df$Cluster) + 1
      working_df$Position_upstream_adj <- working_df$Position_adj+40
      working_df$Position_downstream_adj <- working_df$Position_adj-40
      working_df <- data.frame(Ancestor=anc,Evolved=files[i],working_df)
      
      detail_df <- rbind(detail_df,working_df)
    }
  }
}



write.table(detail_df,'breseq_out.txt',row.names=F)
write.table(summary_df,'breseq_summary.txt',row.names=F)