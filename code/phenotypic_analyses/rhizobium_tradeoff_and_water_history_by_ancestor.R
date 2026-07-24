#here we: characterize tradeoffs in rhizobium partner quality in contemporary wet and dry environments,
#while also examining the impact of rhizobium watering history on partner quality, while controlling for 
#ancestral lineage. Within the Ricks et al manuscript, these analyses correspond to Figure 4 & Tables S10, S11



library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(multcomp)
library(permute)
library(dplyr)

##########################################################################################################
##########################################################################################################
#read in data
##########################################################################################################
##########################################################################################################
#read in plant phenotypic data
data <- read.csv('G:/path/to/greenhouse/phenotypic/data/testphase.csv',header=T)
head(data)

#herbivoroy is coded for the observed herbivory
#plants with no observed herbivory we add a zero to here
data$Herbivory[is.na(data$Herbivory)] <- 0

#create a new variable, root:shoot ratio-ratio between above and belowground biomass
data$RootShoot <- data$Belowground_biomass/data$Aboveground_biomass

#read in mapper file for descendants 
descendant <- read.csv('/path/to/isolation/frequency/data/evolved_strain_identity.csv')

#match phenotypic data to ancestor data
data$ancestor <- NA
for(i in 1:nrow(data)){
  temp <- descendant[descendant$Greenhouse_name==data$Strain[i],]$Ancestor
  if(length(temp)>0){
    data$ancestor[i] <- temp
  }
}

#assign ancestor as a factor
data$ancestor <- as.character(data$ancestor)

##########################################################################################################
##########################################################################################################
#characterize tradeoffs in strain benefits between environments
##########################################################################################################
##########################################################################################################
#we examine the tradeoff in partner quality ofa strain between contemporary wet and contemporary dry environments.
#we focus specifically on the PC axis, our proxy for partner quality. As described elsewhere, our PC is
#composite of our aboveground traits Below, we will additionally include tradeoffs for the other traits. 

#scale variables of interest between 0 and 1
data$Leaf_count_scale <- (data$Leaf_count-min(data$Leaf_count,na.rm=T))/(max(data$Leaf_count,na.rm=T)-min(data$Leaf_count,na.rm=T))
data$Aboveground_biomass_scale <- (data$Aboveground_biomass-min(data$Aboveground_biomass,na.rm=T))/
  (max(data$Aboveground_biomass,na.rm=T)-min(data$Aboveground_biomass,na.rm=T))
data$Height_scale <- (data$Height-min(data$Height,na.rm=T))/(max(data$Height,na.rm=T)-min(data$Height,na.rm=T))

#remove missing values, necessary for PCA
data <- data[!is.na(data$Aboveground_biomass_scale) & !is.na(data$Leaf_count_scale) & !is.na(data$Height_scale),]

#create and extract PCA axes
pc.axis <- (prcomp(cbind(data$Aboveground_biomass_scale,data$Leaf_count_scale,data$Height_scale),scale=T))
data$pc1 <- pc.axis$x[,1]
data$pc2 <- pc.axis$x[,2]

#pc1 is positively correlated with all 3 traits
plot(pc1~Aboveground_biomass_scale,data)
plot(pc1~Leaf_count_scale,data)
plot(pc1~Height_scale,data)


#predict individual strain impacts in both wet and dry predicting individual strain impacts
st_dry.est <- data.frame(emmeans(lm(pc1~Strain,data[data$ContemporaryWater=='Dry',]),~Strain))
st_wet.est <- data.frame(emmeans(lm(pc1~Strain,data[data$ContemporaryWater=='Wet',]),~Strain))

st_dry.est <- tapply(data[data$ContemporaryWater=='Dry',]$pc1,
                     data[data$ContemporaryWater=='Dry',]$Strain,mean)

st_wet.est <- tapply(data[data$ContemporaryWater=='Wet',]$pc1,
                     data[data$ContemporaryWater=='Wet',]$Strain,mean)

merged.strain <- data.frame(Strain=names(st_dry.est),WetStrain=st_wet.est,DryStrain=st_dry.est)

#only looking at data with plant present and no nitrogen fertilizer
#data.sub <- data[data$HistoricPlant=='Present' & data$HistoricNitrogen=='Control',]

#assign historical treatments to each strain
for(i in 1:nrow(merged.strain)){
  temp <- data[data$Strain==merged.strain$Strain[i],][1,]
  merged.strain$HistoricPlant[i] <- temp$HistoricPlant
  merged.strain$HistoricWater[i] <- temp$HistoricWater
  merged.strain$HistoricNitrogen[i] <- temp$HistoricNitrogen
  merged.strain$Ancestor[i] <- temp$ancestor
  
}



#We're only interested in strains that evolved with a plant and under low nitrogen environments
#the other two environments (high nitrogen and plant absent) were associated with lower partner quality
merged.strain <- merged.strain[merged.strain$HistoricPlant=='Present' & merged.strain$HistoricNitrogen=='Control' & !is.na(merged.strain$Ancestor),]

#As we examine these tradeoffs, we need to adjust for the ancestral background. To do this we examine partner quality
#trait variation amongst strains within an ancestral background. While we had 28 original strains background
#in some of these, there were very few replicate isolations. We consequently remove here strains where we have less than 
#2 replicate strains. 
table(merged.strain$Ancestor)#remove samples with only 1 
merged.strain <- merged.strain[!(merged.strain$Ancestor=='262' | merged.strain$Ancestor=='475' |merged.strain$Ancestor=='706'),]


#to adjust for ancestral background, scale partner quality amongst strains from that lineage
st.vec <- unique(merged.strain$Ancestor)
dif.df <- NULL
for(i in 1:length(st.vec)){
  temp <- merged.strain[merged.strain$Ancestor==st.vec[i],]
  temp$DryStrain  <- as.numeric(scale(temp$DryStrain,scale=T))
  temp$WetStrain    <- as.numeric(scale(temp$WetStrain,scale=T))
  dif.df <- rbind(dif.df,temp)
  
}

#correlate partner quality, adjusted by lineage, between wet and dry environments 
plot(DryStrain~WetStrain   ,dif.df)
abline(lm(DryStrain~WetStrain   ,dif.df))
anova(lm(DryStrain~WetStrain,dif.df))
summary(lm(DryStrain~WetStrain,dif.df))




##########################################################################################################
##########################################################################################################
#examine impact of watering history on plant growth, adjusting for lineage
##########################################################################################################
##########################################################################################################
#While the prior models have shown a significant impact on partner quality from the rhizobia's  
#selective environment, these analyses have not accounted for ancestral lineages. Consequently, 
#the observed shifts in traits and partner quality can be due to both in the composition of the 
#ancestral lineages we have (selection on standing variation), as well variation with in 
#these ancestral lineages (likely due to introduced variation). With significant tradeoffs 
#shown above, this suggests significant intra-lineage variation. Below we run models examining 
#partner quality within lineages. These are mixed effects models, using lineage as a random effect




#here we first examine variation within between the historical watering treatments. we use the 
#same set of strains as the above tradeoff analysis, that have sufficient replication
#we look first only at the impact of watering history
#we run these models for every trait, including the PC and a multivariate RDA
subbed.data <- data[data$ancestor %in% st.vec & data$HistoricPlant=='Present' & data$HistoricNitrogen=='Control',]

########################################################################################################
#leaf count
########################################################################################################
anc.leaf_full <- (lmer(Leaf_count~ContemporaryWater*HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anc.leaf_wet <- (lmer(Leaf_count~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Wet' ,]))
anc.leaf_dry <- (lmer(Leaf_count~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Dry' ,]))

anova(anc.leaf_full)
anova(anc.leaf_wet)
anova(anc.leaf_dry)

summary(anc.leaf_full)
summary(anc.leaf_wet)
summary(anc.leaf_dry)

plot(emmeans(anc.leaf_full,~HistoricWater+ContemporaryWater))


########################################################################################################
#height
########################################################################################################
anc.height_full <- (lmer(Height~ContemporaryWater*HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anc.height_wet <- (lmer(Height~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Wet' ,]))
anc.height_dry <- (lmer(Height~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Dry' ,]))

anova(anc.height_full)
anova(anc.height_wet)
anova(anc.height_dry)

summary(anc.height_full)
summary(anc.height_wet)
summary(anc.height_dry)

plot(emmeans(anc.height_full,~HistoricWater+ContemporaryWater))

########################################################################################################
#aboveground biomass
########################################################################################################
anc.above_full <- (lmer(Aboveground_biomass~ContemporaryWater*HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anc.above_wet <- (lmer(Aboveground_biomass~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Wet' ,]))
anc.above_dry <- (lmer(Aboveground_biomass~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Dry' ,]))

anova(anc.above_full)
anova(anc.above_wet)
anova(anc.above_dry)

summary(anc.above_full)
summary(anc.above_wet)
summary(anc.above_dry)

plot(emmeans(anc.above_full,~HistoricWater+ContemporaryWater))

########################################################################################################
#belowground biomass
########################################################################################################
anc.below_full <- (lmer(Belowground_biomass~ContemporaryWater*HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anc.below_wet <- (lmer(Belowground_biomass~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Wet' ,]))
anc.below_dry <- (lmer(Belowground_biomass~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Dry' ,]))

anova(anc.below_full)
anova(anc.below_wet)
anova(anc.below_dry)

summary(anc.below_full)
summary(anc.below_wet)
summary(anc.below_dry)

plot(emmeans(anc.below_full,~HistoricWater+ContemporaryWater))

########################################################################################################
#Nodule count
########################################################################################################
anc.nod_full <- (lmer(Nodule_count~ContemporaryWater*HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anc.nod_wet <- (lmer(Nodule_count~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain) ,subbed.data[subbed.data$ContemporaryWater=='Wet' ,]))
anc.nod_dry <- (lmer(Nodule_count~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Dry' ,]))

anova(anc.nod_full)
anova(anc.nod_wet)
anova(anc.nod_dry)

summary(anc.nod_full)
summary(anc.nod_wet)
summary(anc.nod_dry)

plot(emmeans(anc.nod_full,~HistoricWater+ContemporaryWater))

########################################################################################################
#root-shoot
########################################################################################################
anc.rs_full <- (lmer(RootShoot~ContemporaryWater*HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anc.rs_wet <- (lmer(RootShoot~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Wet' ,]))
anc.rs_dry <- (lmer(RootShoot~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Dry' ,]))

anova(anc.rs_full)
anova(anc.rs_wet)
anova(anc.rs_dry)

summary(anc.rs_full)
summary(anc.rs_wet)
summary(anc.rs_dry)

plot(emmeans(anc.rs_full,~HistoricWater+ContemporaryWater))




########################################################################################################
#pc1 
########################################################################################################
anc.pc_full <- (lmer(pc1~ContemporaryWater*HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anc.pc_wet <- (lmer(pc1~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Wet' ,]))
anc.pc_dry <- (lmer(pc1~HistoricWater + (1|ancestor) + (1|Block) +(1|Strain),subbed.data[subbed.data$ContemporaryWater=='Dry' ,]))

anova(anc.pc_full)
anova(anc.pc_wet)
anova(anc.pc_dry)

summary(anc.pc_full)
summary(anc.pc_wet)
summary(anc.pc_dry)

plot(emmeans(anc.pc_full,~HistoricWater+ContemporaryWater))



########################################################################################################
#multivariate response using rda 
########################################################################################################
#rda models-permutation approach 


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#Analysis for just the contemporary wet
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
data.sub.wet <- subbed.data[subbed.data$ContemporaryWater=='Wet',]


#create matrix of traits of interest
traits_matrix <- data.sub.wet[,colnames(data.sub.wet) %in% c('Leaf_count_scale','Height_scale','Aboveground_biomass_scale')]

#calculate the observed f-statistic from the RDA
obs_rda <- rda(traits_matrix ~ HistoricWater+ Condition(ancestor) + Condition(Block) + Condition(Herbivory),data=data.sub.wet)

anova(obs_rda, by = "term")
obs_anova <- anova(obs_rda, by = "term", permutations = 0)
#remove residual row
term_names <- rownames(obs_anova)[1:(nrow(obs_anova)-1)] 
obs_Fs <- obs_anova$F[1:length(term_names)] 
names(obs_Fs) <- term_names

#map strains to their treatments
strain_map <- data.sub.wet %>% 
  select(Strain, HistoricWater) %>% 
  distinct()

#set up the custom permutation
n_perms <- 999
#create a matrix to store the null Fs (Rows = permutations, Columns = model terms)
null_F_matrix <- matrix(NA, nrow = n_perms, ncol = length(obs_Fs))
colnames(null_F_matrix) <- term_names

for(i in 1:n_perms) {
  #shuffle the strain IDs, leaving the treatment combinations intact
  #assigns a whole strain to a random evolutionary history block
  shuffled_map <- strain_map
  shuffled_map$Strain <- sample(shuffled_map$Strain)
  
  #map the new fake histories back to the unbalanced replicates
  df_shuffled <- data.sub.wet %>%
    select(-HistoricWater) %>% # Remove the true treatments
    left_join(shuffled_map, by = "Strain") # Attach the shuffled ones
  
  #run the RDA on the randomized data
  perm_rda <- rda(traits_matrix ~ HistoricWater + Condition(ancestor) + Condition(Block) + Condition(Herbivory), data = df_shuffled)
  
  #extract the null F-statistics for all terms
  perm_anova <- anova(perm_rda, by = "term", permutations = 0)
  null_F_matrix[i, ] <- perm_anova$F[1:length(obs_Fs)]
  print(i)
}

#calculate custom p-values for every term
p_values_wet <- numeric(length(obs_Fs))
names(p_values_wet) <- term_names

for(j in 1:length(obs_Fs)) {
  #for each term, how many null Fs are greater than or equal to the observed
  p_values_wet[j] <- sum(null_F_matrix[, j] >= obs_Fs[j]) / (n_perms + 1)
}
p_values_wet



#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#Analysis for just the contemporary dry
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
data.sub.dry <- subbed.data[subbed.data$ContemporaryWater=='Dry',]


#create matrix of traits of interest
traits_matrix <- data.sub.dry[,colnames(data.sub.dry) %in% c('Leaf_count_scale','Height_scale','Aboveground_biomass_scale')]

#calculate the observed f-statistic from the RDA
obs_rda <- rda(traits_matrix ~ HistoricWater+ Condition(ancestor) + Condition(Block) + Condition(Herbivory),data=data.sub.dry)
anova(obs_rda, by = "term")
obs_anova <- anova(obs_rda, by = "term", permutations = 0)
#remove residual row
term_names <- rownames(obs_anova)[1:(nrow(obs_anova)-1)] 
obs_Fs <- obs_anova$F[1:length(term_names)] 
names(obs_Fs) <- term_names

#map strains to their treatments
strain_map <- data.sub.dry %>% 
  select(Strain, HistoricWater) %>% 
  distinct()

#set up the custom permutation
n_perms <- 999
#create a matrix to store the null Fs (Rows = permutations, Columns = model terms)
null_F_matrix <- matrix(NA, nrow = n_perms, ncol = length(obs_Fs))
colnames(null_F_matrix) <- term_names

for(i in 1:n_perms) {
  #shuffle the strain IDs, leaving the treatment combinations intact
  #assigns a whole strain to a random evolutionary history block
  shuffled_map <- strain_map
  shuffled_map$Strain <- sample(shuffled_map$Strain)
  
  #map the new fake histories back to the unbalanced replicates
  df_shuffled <- data.sub.dry %>%
    select(-HistoricWater) %>% # Remove the true treatments
    left_join(shuffled_map, by = "Strain") # Attach the shuffled ones
  
  #run the RDA on the randomized data
  perm_rda <- rda(traits_matrix ~ HistoricWater + Condition(ancestor) + Condition(Block) + Condition(Herbivory), data = df_shuffled)
  
  #extract the null F-statistics for all terms
  perm_anova <- anova(perm_rda, by = "term", permutations = 0)
  null_F_matrix[i, ] <- perm_anova$F[1:length(obs_Fs)]
  print(i)
}

#calculate custom p-values for every term
p_values_dry <- numeric(length(obs_Fs))
names(p_values_wet) <- term_names

for(j in 1:length(obs_Fs)) {
  #for each term, how many null Fs are greater than or equal to the observed
  p_values_dry[j] <- sum(null_F_matrix[, j] >= obs_Fs[j]) / (n_perms + 1)
}
p_values_dry


#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#Analysis for across both watering environments
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#create matrix of traits of interest
traits_matrix <- subbed.data[,colnames(subbed.data) %in% c('Leaf_count_scale','Height_scale','Aboveground_biomass_scale')]

#calculate the observed f-statistic from the RDA
obs_rda <- rda(traits_matrix ~ ContemporaryWater*HistoricWater+ Condition(ancestor) + Condition(Herbivory),data=subbed.data)
anova(obs_rda, by = "term")
obs_anova <- anova(obs_rda, by = "term", permutations = 0)
#remove residual row
term_names <- rownames(obs_anova)[1:(nrow(obs_anova)-1)] 
obs_Fs <- obs_anova$F[1:length(term_names)] 
names(obs_Fs) <- term_names

#map strains to their treatments
strain_map <- subbed.data %>% 
  select(Strain, HistoricWater) %>% 
  distinct()

#set up the custom permutation
n_perms <- 999
#create a matrix to store the null Fs (Rows = permutations, Columns = model terms)
null_F_matrix <- matrix(NA, nrow = n_perms, ncol = length(obs_Fs))
colnames(null_F_matrix) <- term_names

for(i in 1:n_perms) {
  #shuffle the strain IDs, leaving the treatment combinations intact
  #assigns a whole strain to a random evolutionary history block
  shuffled_map <- strain_map
  shuffled_map$Strain <- sample(shuffled_map$Strain)
  
  #map the new fake histories back to the unbalanced replicates
  df_shuffled <- subbed.data %>%
    select(-HistoricWater) %>% # Remove the true treatments
    left_join(shuffled_map, by = "Strain") # Attach the shuffled ones
  
  #run the RDA on the randomized data
  perm_rda <- rda(traits_matrix ~ ContemporaryWater*HistoricWater + Condition(ancestor) +  Condition(Herbivory), data = df_shuffled)
  
  #extract the null F-statistics for all terms
  perm_anova <- anova(perm_rda, by = "term", permutations = 0)
  null_F_matrix[i, ] <- perm_anova$F[1:length(obs_Fs)]
  print(i)
}

#calculate custom p-values for every term
p_values_mixed <- numeric(length(obs_Fs))
names(p_values_mixed) <- term_names

for(j in 1:length(obs_Fs)) {
  #for each term, how many null Fs are greater than or equal to the observed
  p_values_mixed[j] <- sum(null_F_matrix[, j] >= obs_Fs[j]) / (n_perms + 1)
}
p_values_mixed

