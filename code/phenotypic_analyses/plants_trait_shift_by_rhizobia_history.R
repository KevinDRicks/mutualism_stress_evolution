#here we examine the impact of  water, nitrogen, and plant history on rhizobium partner quality
#we examine multiple plant traits of interest, as well as using multivariate approaches (pc for 
#variable compression and rda for full multivariate analysis). Within the Ricks et al manuscript, 
#these analyses correspond to Figure 1, S2, S3 & Tables S1, S2


library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(multcomp)
library(permute)
library(dplyr)

##########################################################################################################
#read in and process phenotype data from test phase of experiment
#(greenhouse experiment where we evaluated phenotypic effects of all evolved strains, under both wet and dry environments)
##########################################################################################################
data <- read.csv('G:/path/to/greenhouse/phenotypic/data/testphase.csv',header=T)
head(data)

#herbivoroy is coded for the observed herbivory
#plants with no observed herbivory we add a zero to here
data$Herbivory[is.na(data$Herbivory)] <- 0

#create a new variable, root:shoot ratio-ratio between above and belowground biomass
data$RootShoot <- data$Belowground_biomass/data$Aboveground_biomass

########################################################################################################
########################################################################################################
#univariate analysis of plant traits
########################################################################################################
########################################################################################################

#linear mixed effects models evaluate impacts of contemporary water and interactions with strain 
#evolution historyin impacting individual plant traits. We include Herbivory as a covariate and
#greenhouse block and the strain as random effects. For every trait, we run a full model, 
#across both contemporary watering treatment. However, as interpreting 4-way interactions can be
#challenging (ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant) we additionally run
#models for each contemporary watering environment
#Quick note: some models produced singular fits, as the strain genotype term explained little to no 
#additional variance. While this random effect could be dropped, we retained it in all models to 
#correctly reflect the study design, as replicate plants within a genotype are not statistically 
#independent


########################################################################################################
#Leaf count
########################################################################################################
leaf_full <- lmer(Leaf_count~ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
           data=data)
leaf_wet <- lmer(Leaf_count~HistoricWater*HistoricPlant*HistoricNitrogen + Herbivory + (1|Block) + (1|Strain), 
                  data=data[data$ContemporaryWater=='Wet',])
leaf_dry <- lmer(Leaf_count~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                 data=data[data$ContemporaryWater=='Dry',])
anova(leaf_full)
anova(leaf_wet)
anova(leaf_dry)
summary(leaf_full)
summary(leaf_wet)
summary(leaf_dry)

plot(emmeans(leaf_wet,~HistoricWater*HistoricNitrogen*HistoricPlant))
plot(emmeans(leaf_dry,~HistoricWater*HistoricNitrogen*HistoricPlant))


########################################################################################################
#Height
########################################################################################################
height_full <- lmer(Height~ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                  data=data)
height_wet <- lmer(Height~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                 data=data[data$ContemporaryWater=='Wet',])
height_dry <- lmer(Height~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                 data=data[data$ContemporaryWater=='Dry',])
anova(height_full)
anova(height_wet)
anova(height_dry)
summary(height_full)
summary(height_wet)
summary(height_dry)

plot(emmeans(height_wet,~HistoricWater*HistoricNitrogen*HistoricPlant))
plot(emmeans(height_dry,~HistoricWater*HistoricNitrogen*HistoricPlant))


########################################################################################################
#Aboveground biomass
########################################################################################################
above_full <- lmer(Aboveground_biomass~ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                    data=data)
above_wet <- lmer(Aboveground_biomass~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                   data=data[data$ContemporaryWater=='Wet',])
above_dry <- lmer(Aboveground_biomass~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                   data=data[data$ContemporaryWater=='Dry',])
anova(above_full)
anova(above_wet)
anova(above_dry)
summary(above_full)
summary(above_wet)
summary(above_dry)

plot(emmeans(above_wet,~HistoricWater*HistoricNitrogen*HistoricPlant))
plot(emmeans(above_dry,~HistoricWater*HistoricNitrogen*HistoricPlant))


########################################################################################################
#Belowground biomass
########################################################################################################
below_full <- lmer(Belowground_biomass~ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                    data=data)
below_wet <- lmer(Belowground_biomass~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                   data=data[data$ContemporaryWater=='Wet',])
below_dry <- lmer(Belowground_biomass~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                   data=data[data$ContemporaryWater=='Dry',])
anova(below_full)
anova(below_wet)
anova(below_dry)
summary(below_full)
summary(below_wet)
summary(below_dry)

plot(emmeans(below_wet,~HistoricWater*HistoricNitrogen*HistoricPlant))
plot(emmeans(below_dry,~HistoricWater*HistoricNitrogen*HistoricPlant))


########################################################################################################
#rootshoot ratio
########################################################################################################
rs_full <- lmer(RootShoot~ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                   data=data)
rs_wet <- lmer(RootShoot~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                  data=data[data$ContemporaryWater=='Wet',])
rs_dry <- lmer(RootShoot~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                  data=data[data$ContemporaryWater=='Dry',])
anova(rs_full)
anova(rs_wet)
anova(rs_dry)
summary(rs_full)
summary(rs_wet)
summary(rs_dry)

plot(emmeans(below_wet,~HistoricWater*HistoricNitrogen*HistoricPlant))
plot(emmeans(below_dry,~HistoricWater*HistoricNitrogen*HistoricPlant))
########################################################################################################
#Nodule count
########################################################################################################
nodule_full <- lmer(Nodule_count~ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain) + (1|NoduleCounter), 
                    data=data)
nodule_wet <- lmer(Nodule_count~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain) + (1|NoduleCounter), 
                   data=data[data$ContemporaryWater=='Wet',])
nodule_dry <- lmer(Nodule_count~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain) , 
                   data=data[data$ContemporaryWater=='Dry',])
anova(nodule_full)
anova(nodule_wet)
anova(nodule_dry)
summary(nodule_full)
summary(nodule_wet)
summary(nodule_dry)

plot(emmeans(nodule_wet,~HistoricWater*HistoricNitrogen*HistoricPlant))
plot(emmeans(nodule_dry,~HistoricWater*HistoricNitrogen*HistoricPlant))


########################################################################################################
########################################################################################################
#multivariate analysis of plant traits
########################################################################################################
########################################################################################################
#we also evaluate the multivariate suite of traits in unison. We specfically examine the three
#aboveground traits, including aboveground biomass, leaf count, and height. Together these
#are all associated with growth rate, and our proxy for rhizobia partner quality. To this end, 
#we use two approaches: 1) linear models on pc axes and 2) explicitly multivariate models in 
#the form of an rda. 

########################################################################################################
#analysis of pc axes
########################################################################################################
#With ~2000 plants and 5 measured plant traits, there were some failures in coordination and failured germination.
#consequently, a small number of samples have missing data. for pca, there can't be missing data
#we remove sapmles that are missing one of the 3 traits of interest
data.sub <- data


#scale variables of interest between 0 and 1
data.sub$Leaf_count <- (data.sub$Leaf_count-min(data.sub$Leaf_count,na.rm=T))/(max(data.sub$Leaf_count,na.rm=T)-min(data.sub$Leaf_count,na.rm=T))
data.sub$Aboveground_biomass <- (data.sub$Aboveground_biomass-min(data.sub$Aboveground_biomass,na.rm=T))/
  (max(data.sub$Aboveground_biomass,na.rm=T)-min(data.sub$Aboveground_biomass,na.rm=T))
data.sub$Height <- (data.sub$Height-min(data$Height,na.rm=T))/(max(data.sub$Height,na.rm=T)-min(data$Height,na.rm=T))

data.sub <- data.sub[!is.na(data.sub$Aboveground_biomass) & !is.na(data.sub$Leaf_count) & !is.na(data.sub$Height),]


#genereate pc axes
pc.axis <- (prcomp(cbind(data.sub$Aboveground_biomass,data.sub$Leaf_count,data.sub$Height),scale=T))

#variance explained by pc axes
summary(pc.axis)
summary(pc.axis)$importance[2,]

data.sub$pc1 <- pc.axis$x[,1]
data.sub$pc2 <- pc.axis$x[,2]
plot(pc1~Aboveground_biomass,data.sub)
plot(pc1~Height,data.sub)
plot(pc1~Leaf_count,data.sub)

plot(pc2~Aboveground_biomass,data.sub)
plot(pc2~Height,data.sub)
plot(pc2~Leaf_count,data.sub)
#pc1 captures the majority of the variation, and is positively correlated with all 3 variables
#we use this as our proxy for partner quality moving forward



#run models
pc1.full <- lmer(pc1~ContemporaryWater*HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
           data=data.sub)
pc1.wet <- lmer(pc1~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                 data=data.sub[data.sub$ContemporaryWater=='Wet',])
pc1.dry <- lmer(pc1~HistoricNitrogen*HistoricWater*HistoricPlant + Herbivory + (1|Block) + (1|Strain), 
                data=data.sub[data.sub$ContemporaryWater=='Dry',])
anova(pc1.full)
anova(pc1.wet)
anova(pc1.dry)

summary(pc1.full)
summary(pc1.wet)
summary(pc1.dry)

plot(emmeans(pc1.wet,~HistoricWater*HistoricNitrogen*HistoricPlant))
plot(emmeans(pc1.dry,~HistoricWater*HistoricNitrogen*HistoricPlant))



########################################################################################################
#multivariate rda analysis
########################################################################################################
#running the rda, we need to account for pseudo replication within each evolved strain. Unlike univariate 
#mixed models above, we've can't put these as a random effect. rda can generate stats on the inputted 
#variables by permuting data. standard approaches for correcting for pseudo replicate within rda models 
#use controlled permutation. In this case, samples within a single strain are moved together in the 
#permutation. We run a custom permutation approach below 




#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#Analysis for just the contemporary wet
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
data.sub.wet <- data.sub[data.sub$ContemporaryWater=='Wet',]


#create matrix of traits of interest
traits_matrix <- data.sub.wet[,colnames(data.sub.wet) %in% c('Leaf_count','Height','Aboveground_biomass')]

#calculate the observed f-statistic from the RDA
obs_rda <- rda(traits_matrix ~ HistoricWater*HistoricNitrogen*HistoricPlant+ Condition(Block) + Condition(Herbivory),data=data.sub.wet)

anova(obs_rda, by = "term")
obs_anova <- anova(obs_rda, by = "term", permutations = 0)
#remove residual row
term_names <- rownames(obs_anova)[1:(nrow(obs_anova)-1)] 
obs_Fs <- obs_anova$F[1:length(term_names)] 
names(obs_Fs) <- term_names

#map strains to their treatments
strain_map <- data.sub.wet %>% 
  select(Strain, HistoricWater, HistoricNitrogen, HistoricPlant) %>% 
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
    select(-HistoricWater, -HistoricNitrogen, -HistoricPlant) %>% # Remove the true treatments
    left_join(shuffled_map, by = "Strain") # Attach the shuffled ones
  
  #run the RDA on the randomized data
  perm_rda <- rda(traits_matrix ~ HistoricWater * HistoricNitrogen * HistoricPlant + Condition(Block) + Condition(Herbivory), data = df_shuffled)
  
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
data.sub.dry <- data.sub[data.sub$ContemporaryWater=='Dry',]


#create matrix of traits of interest
traits_matrix <- data.sub.dry[,colnames(data.sub.dry) %in% c('Leaf_count','Height','Aboveground_biomass')]

#calculate the observed f-statistic from the RDA
obs_rda <- rda(traits_matrix ~ HistoricWater*HistoricNitrogen*HistoricPlant+ Condition(Block) + Condition(Herbivory),data=data.sub.dry)
anova(obs_rda, by = "term")
obs_anova <- anova(obs_rda, by = "term", permutations = 0)
#remove residual row
term_names <- rownames(obs_anova)[1:(nrow(obs_anova)-1)] 
obs_Fs <- obs_anova$F[1:length(term_names)] 
names(obs_Fs) <- term_names

#map strains to their treatments
strain_map <- data.sub.dry %>% 
  select(Strain, HistoricWater, HistoricNitrogen, HistoricPlant) %>% 
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
    select(-HistoricWater, -HistoricNitrogen, -HistoricPlant) %>% # Remove the true treatments
    left_join(shuffled_map, by = "Strain") # Attach the shuffled ones
  
  #run the RDA on the randomized data
  perm_rda <- rda(traits_matrix ~ HistoricWater * HistoricNitrogen * HistoricPlant + Condition(Block) + Condition(Herbivory), data = df_shuffled)
  
  #extract the null F-statistics for all terms
  perm_anova <- anova(perm_rda, by = "term", permutations = 0)
  null_F_matrix[i, ] <- perm_anova$F[1:length(obs_Fs)]
  print(i)
}

#calculate custom p-values for every term
p_values_dry <- numeric(length(obs_Fs))
names(p_values_dry) <- term_names

for(j in 1:length(obs_Fs)) {
  #for each term, how many null Fs are greater than or equal to the observed
  p_values_dry[j] <- sum(null_F_matrix[, j] >= obs_Fs[j]) / (n_perms + 1)
}
p_values_dry




#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
#Analysis for across both watering environments
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

#create matrix of traits of interest
traits_matrix <- data.sub[,colnames(data.sub) %in% c('Leaf_count','Height','Aboveground_biomass')]

#calculate the observed f-statistic from the RDA
obs_rda <- rda(traits_matrix ~ ContemporaryWater*HistoricWater*HistoricNitrogen*HistoricPlant + Condition(Herbivory),data=data.sub)

anova(obs_rda, by = "term")
obs_anova <- anova(obs_rda, by = "term", permutations = 0)
#remove residual row
term_names <- rownames(obs_anova)[1:(nrow(obs_anova)-1)] 
obs_Fs <- obs_anova$F[1:length(term_names)] 
names(obs_Fs) <- term_names

#map strains to their treatments
strain_map <- data.sub %>% 
  select(Strain, HistoricWater, HistoricNitrogen, HistoricPlant) %>% 
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
  df_shuffled <- data.sub %>%
    select(-HistoricWater, -HistoricNitrogen, -HistoricPlant) %>% # Remove the true treatments
    left_join(shuffled_map, by = "Strain") # Attach the shuffled ones
  
  #run the RDA on the randomized data
  perm_rda <- rda(traits_matrix ~ ContemporaryWater * HistoricWater * HistoricNitrogen * HistoricPlant + Condition(Herbivory), data = df_shuffled)
  
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



