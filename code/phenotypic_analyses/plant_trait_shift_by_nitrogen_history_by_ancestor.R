#here we examine the impact of rhizobium nitrogen history on partner quality, while controlling for 
#ancestral lineage. Within the Ricks et al manuscript, these analyses correspond to Figure S12 & Tables S12


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
#examine impact of nitrogen history on plant growth, adjusting for lineage
##########################################################################################################
##########################################################################################################
#Below we run models examining partner quality within lineages. These are mixed effects models, using 
#lineage as a random effect

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

#strongest effects for nitrogen were found when evolved with a plant and in wet environments, and grown in contemporary wet
#subset to these group
subbed.data <- data[data$ContemporaryWater=='Wet' & 
                      data$HistoricWater=='Wet' & data$HistoricPlant=='Present' & !is.na(data$ancestor),]

#subset to only include ancestral lineages with replicates in both control and nitrogen groups
control.rep <- data.frame(subbed.data[subbed.data$HistoricNitrogen=='Control',] %>%
  group_by(ancestor) %>%  summarise(n_strains = n_distinct(Strain)))
nitrogen.rep <- data.frame(subbed.data[subbed.data$HistoricNitrogen=='Nitrogen',] %>%
                            group_by(ancestor) %>%  summarise(n_strains = n_distinct(Strain)))

subbed.data <- subbed.data[subbed.data$ancestor %in% intersect(control.rep$ancestor,nitrogen.rep$ancestor),]

#here we first examine variation within between the historical nitrogen treatments.
#we run these models for every trait, including the PC and a multivariate RDA
########################################################################################################
#leaf count
########################################################################################################
anc.leaf <- (lmer(Leaf_count~HistoricNitrogen + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anova(anc.leaf)
summary(anc.leaf)
plot(emmeans(anc.leaf,~HistoricNitrogen))


########################################################################################################
#height
########################################################################################################
anc.height <- (lmer(Height~HistoricNitrogen + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anova(anc.height)
summary(anc.height)
plot(emmeans(anc.height,~HistoricNitrogen))


########################################################################################################
#aboveground biomass
########################################################################################################
anc.above <- (lmer(Aboveground_biomass~HistoricNitrogen + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anova(anc.above)
summary(anc.above)
plot(emmeans(anc.above,~HistoricNitrogen))


########################################################################################################
#belowground biomass
########################################################################################################
anc.below <- (lmer(Belowground_biomass~HistoricNitrogen + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anova(anc.below)
summary(anc.below)
plot(emmeans(anc.below,~HistoricNitrogen))


########################################################################################################
#Nodule count
########################################################################################################
anc.nod <- (lmer(Nodule_count~HistoricNitrogen + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anova(anc.nod)
summary(anc.nod)
plot(emmeans(anc.nod,~HistoricNitrogen))


########################################################################################################
#root-shoot
########################################################################################################
anc.rs <- (lmer(RootShoot~HistoricNitrogen + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anova(anc.rs)
summary(anc.rs)
plot(emmeans(anc.rs,~HistoricNitrogen))





########################################################################################################
#pc1 
########################################################################################################
anc.pc1 <- (lmer(pc1~HistoricNitrogen + (1|ancestor) + (1|Block) +(1|Strain),subbed.data))
anova(anc.pc1)
summary(anc.pc1)
plot(emmeans(anc.pc1,~HistoricNitrogen))



########################################################################################################
#multivariate response using rda 
########################################################################################################
#rda models-permutation approach 

#create matrix of traits of interest
traits_matrix <- subbed.data[,colnames(subbed.data) %in% c('Leaf_count_scale','Height_scale','Aboveground_biomass_scale')]

#calculate the observed f-statistic from the RDA
obs_rda <- rda(traits_matrix ~ HistoricNitrogen+ Condition(ancestor) + Condition(Block) + Condition(Herbivory),data=subbed.data)
anova(obs_rda, by = "term")
obs_anova <- anova(obs_rda, by = "term", permutations = 0)
#remove residual row
term_names <- rownames(obs_anova)[1:(nrow(obs_anova)-1)] 
obs_Fs <- obs_anova$F[1:length(term_names)] 
names(obs_Fs) <- term_names

#map strains to their treatments
strain_map <- subbed.data %>% 
  select(Strain, HistoricNitrogen) %>% 
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
    select(-HistoricNitrogen) %>% # Remove the true treatments
    left_join(shuffled_map, by = "Strain") # Attach the shuffled ones
  
  #run the RDA on the randomized data
  perm_rda <- rda(traits_matrix ~ HistoricNitrogen+ Condition(ancestor) + Condition(Block) +  Condition(Herbivory), data = df_shuffled)
  
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

