#analysis of gene presence absence variation among treatments
#we do two pairwise comparisons. First between the wet and dry, low nitrogen treatments (pwc vs pdc). 
#secondly, between the low and high nitrogen treatments, from the wet treatments (pwc vs pwn).
#we constrain this nitrogen treatments to solely the wet environments as the largest effects of nitrogen
#could be found in these treatments. Within the Ricks et al manuscript, 
#these analyses correspond to Figure S8, and Tables S5 and S6


########################################################################################################
#comparison between wet and dry  (pwc vs pdc)
########################################################################################################
#read in data
gpi <- read.csv('path/to/presence/absence/gene/variation/data/pav_pwc_vs_pdc.csv')
head(gpi)
colnames(gpi)
dim(gpi)
#remove genes that found in every isolate
#these are the background/core genes
#there are 38 PDC and 40 PWC isolates
gpi <- gpi[!(gpi$PDC_freq==38 & gpi$PWC_freq==40),]
dim(gpi)


#loop through a fishers test across all genes 
fp.vec <- NULL
odds.vec <- NULL
for(i in 1:nrow(gpi)){
  pdc.pres <- gpi$PDC_freq[i]
  pwc.pres <- gpi$PWC_freq[i]
  pdc.abs <- 38-pdc.pres
  pwc.abs <- 40-pwc.pres
  M <- as.table(rbind(c(pdc.pres,pdc.abs),c(pwc.pres,pwc.abs)))
  test <- fisher.test(M)
  fp.vec <- c(fp.vec,test$p.value)
  odds.vec <- c(odds.vec,test$estimate)
  
}

#run fdr adjustment
padj <- p.adjust(fp.vec,method='fdr')
sum(padj<0.05)
gpi$pval <- padj

gpi[padj<0.05,]
padj[padj<0.05]
table(gpi[padj<0.05,]$Replicon)


########################################################################################################
#comparison between low and high nitrogen  (pwc vs pwn)
########################################################################################################

#read in data
gpi <- read.csv('path/to/presence/absence/gene/variation/data/pav_pwc_vs_pwn.csv')
head(gpi)
colnames(gpi)
dim(gpi)
#remove genes that found in every isolate
#these are the background/core genes
#there are 40 PWC and 38 PWN isolates
gpi <- gpi[!(gpi$PWC_freq==40 & gpi$PWN_freq==38),]
dim(gpi)


fp.vec <- NULL
odds.vec <- NULL
for(i in 1:nrow(gpi)){
  pwn.pres <- gpi$PWN_freq[i]
  pwc.pres <- gpi$PWC_freq[i]
  pwn.abs <- max(gpi$PWN_freq)-pwn.pres
  pwc.abs <- max(gpi$PWC_freq)-pwc.pres
  M <- as.table(rbind(c(pwn.pres,pwn.abs),c(pwc.pres,pwc.abs)))
  test <- fisher.test(M)
  fp.vec <- c(fp.vec,test$p.value)
  odds.vec <- c(odds.vec,test$estimate)
}



#run fdr adjustment
padj <- p.adjust(fp.vec,method='fdr')
sum(padj<0.05)
gpi$pval <- padj

gpi[padj<0.05,]
padj[padj<0.05]
table(gpi[padj<0.05,]$Replicon)

