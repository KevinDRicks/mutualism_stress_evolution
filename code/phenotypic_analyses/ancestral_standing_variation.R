#here we examine the variation in the isolation frequence of our ancestral rhizobia strains, by treatment
#we then evaluate the partner quality in our ancestral rhizobia strains, in both
#contemporary wet and contemporary dry environments
#and relate this ancestral partner quality to their isolation frequency 
#Within the Ricks et al manuscript, these analyses correspond to Figures 2A-C, S1, & S6 &
#Table S3 & S4

library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(ciTools)
library(car)
library(MuMIn)

###########################################################################################
###########################################################################################
#read in data
###########################################################################################
###########################################################################################
#read in ancestor isolation frequence

recovered_strains <- read.csv('/path/to/isolation/frequency/data/evolved_strain_identity.csv')
anc <- unique(recovered_strains$Ancestor)

#create dataframe for breakdown of strains by treatment
anc.d <- NULL
for(i in 1:length(anc)){
  pdc <- nrow(recovered_strains[recovered_strains$Ancestor==anc[i] & recovered_strains$TreatmentCode=='PDC',])
  pdn <- nrow(recovered_strains[recovered_strains$Ancestor==anc[i] & recovered_strains$TreatmentCode=='PDN',])
  pwc <- nrow(recovered_strains[recovered_strains$Ancestor==anc[i] & recovered_strains$TreatmentCode=='PWC',])
  pwn <- nrow(recovered_strains[recovered_strains$Ancestor==anc[i] & recovered_strains$TreatmentCode=='PWN',])
  anc.d <- rbind(anc.d,c(pdc,pdn,pwc,pwn))
}

lab.desc <- data.frame(Strain=anc,anc.d)
colnames(lab.desc)[c(2:5)] <- c('PDC','PDN','PWC','PWN')
lab.desc


#run chisq looking at difference between treatments
#here we look at difference between PDC (plant, drought, control) vs PWC (plant, wet, control)
contingency_table <- as.matrix(lab.desc[, c(2,4)])  # Select only PDC and PWC columns
rownames(contingency_table) <- rownames(lab.desc)   # Set row names as strain names
# Remove rows where both PDC and PWC are zero
contingency_table <- contingency_table[rowSums(contingency_table) > 0, ]
chisq.test(contingency_table, simulate.p.value = TRUE,B=100000)



#reads in ancestor phenotype data
anc.data <- read.csv('/path/to/ancestral/plant/phenotype/ancestral_quality.csv')


#remove plants from sterile treatments, those that died, and those with missing data
anc.data <- anc.data[!(anc.data$Strain=='Sterile' | anc.data$Notes=='Missing roots' | anc.data$Notes=='Dead'),]
anc.data <- anc.data[!is.na(anc.data$Aboveground_biomass) & !is.na(anc.data$Leaf_count) & !is.na(anc.data$Height),]


#scale variables of interest between 0 and 1
anc.data$Leaf_count <- (anc.data$Leaf_count-min(anc.data$Leaf_count,na.rm=T))/(max(anc.data$Leaf_count,na.rm=T)-min(anc.data$Leaf_count,na.rm=T))
anc.data$Aboveground_biomass <- (anc.data$Aboveground_biomass-min(anc.data$Aboveground_biomass,na.rm=T))/(max(anc.data$Aboveground_biomass,na.rm=T)-min(anc.data$Aboveground_biomass,na.rm=T))
anc.data$Height <- (anc.data$Height-min(anc.data$Height,na.rm=T))/(max(anc.data$Height,na.rm=T)-min(anc.data$Height,na.rm=T))

#create and extract PCA axes
pc.axis <- (prcomp(cbind(anc.data$Aboveground_biomass,anc.data$Leaf_count,anc.data$Height),scale=T))
anc.data$pc1 <- pc.axis$x[,1]
anc.data$pc2 <- pc.axis$x[,2]

#pc1 is positively correlated with all 3 traits
plot(pc1~Aboveground_biomass,anc.data)
plot(pc1~Leaf_count,anc.data)
plot(pc1~Height,anc.data)


#estimate strain effects on composite pc axis
wet.mod <- lmer((pc1)~Strain + (1|Tray),anc.data[anc.data$Water=='Wet',])
est.wet <- data.frame(emmeans(wet.mod,~Strain))
dry.mod <- lmer((pc1 )~Strain + (1|Tray),anc.data[anc.data$Water=='Dry',])
est.dry <- data.frame(emmeans(dry.mod,~Strain))

anova(wet.mod)
anova(dry.mod)

plot(emmeans(wet.mod,~Strain))
plot(emmeans(dry.mod,~Strain))


#we want to correlate strain isolation frequencies in the pwc and pdc treatments (plant control, in either wet or dry)
#with their strain benefits from this ancestor experiment


#append the strain frequencies to our dataframes with partner qualitly
est.wet$PWC <-0
for(i in 1:nrow(est.wet)){
  temp.work <- lab.desc[lab.desc$Strain==est.wet$Strain[i],]
  if(nrow(temp.work)!=0){
    est.wet$PWC[i] <- temp.work$PWC
  }
}

est.dry$PDC <-0
for(i in 1:nrow(est.dry)){
  temp.work <- lab.desc[lab.desc$Strain==est.dry$Strain[i],]
  if(nrow(temp.work)!=0){
    est.dry$PDC[i] <- temp.work$PDC
  }
}



###############################################################################################################
#correlate partner quality in wet with isolation frequency in wet
###############################################################################################################
#base model
mod1 <- (glm(PWC~emmean,est.wet,family='poisson'))
anova(mod1)

#create prediction interval
x.predict <- seq(-5,max(est.wet$emmean)*1.1,.01)
add_ci(est.wet, mod1, names = c("lcb", "ucb"), alpha = 0.05)
pred.int <- add_ci(data.frame(emmean=x.predict), mod1, names = c("lcb", "ucb"), alpha = 0.05)


plot(PWC~emmean,est.wet,ylim=c(0,7))
lines(pred~emmean,pred.int)
lines(lcb~emmean,pred.int)
lines(ucb~emmean,pred.int)


#compare models of PWC frequency using partner quality in wet vs partner quality in dr
wet.mod <- glm(est.wet$PWC~est.wet$emmean,family='poisson')
dry.mod <- glm(est.wet$PWC~est.dry$emmean,family='poisson')
anova(wet.mod)
anova(dry.mod)

#generally see slightly more variance explained in wet isolation frequency by partner quality in wet
1-(wet.mod$deviance/wet.mod$null.deviance)
1-(dry.mod$deviance/dry.mod$null.deviance)


###############################################################################################################
#correlate partner quality in dry with isolation frequency in dry
###############################################################################################################
#base model
mod1 <- (glm(PDC~emmean,est.dry,family='poisson'))
anova(mod1)

#create prediction interval
x.predict <- seq(-5,0,.01)
add_ci(est.wet, mod1, names = c("lcb", "ucb"), alpha = 0.05)
pred.int <- add_ci(data.frame(emmean=x.predict), mod1, names = c("lcb", "ucb"), alpha = 0.05)


plot(PDC~emmean,est.dry,ylim=c(0,7))
lines(pred~emmean,pred.int)
lines(lcb~emmean,pred.int)
lines(ucb~emmean,pred.int)


#compare models of PDC frequency using partner quality in wet vs partner quality in dr
wet.mod <- glm(est.dry$PDC~est.wet$emmean,family='poisson')
dry.mod <- glm(est.dry$PDC~est.dry$emmean,family='poisson')
anova(wet.mod)
anova(dry.mod)

#generally see slightly more variance explained in dry isolation frequency by partner quality in dry
1-(wet.mod$deviance/wet.mod$null.deviance)
1-(dry.mod$deviance/dry.mod$null.deviance)



###############################################################################################################
#model comparison
###############################################################################################################
#we evaluate if partner quality in the dry versus the wet environment better explain isolation frequency

#rescale partner quality within each environment, so can be comparable
est.dry$emmean <- as.numeric(scale(est.dry$emmean,scale=T))
est.wet$emmean <- as.numeric(scale(est.wet$emmean,scale=T))

#create a new dataframe merging data
#this includes strain identity (strain column), its isolation frequency in either wet or dry (IsoFreq and IsoEnv columns)
#and and the partner quality in wet or dry environments. We however code partner quality as either matched or mismatched
#partner quality in the dry environment is considered a matched environment for isolation frequencies from the dry
#with this approach we can compare matched vs mismatched partner quality as predictors of isolation frequency 
mixed.df <- rbind(data.frame(Strain=est.dry$Strain,IsoFreq=est.dry$PDC,IsoEnv='Dry',
                             Matched_PQ=est.dry$emmean,MisMatched_PQ=est.wet$emmean),
                  data.frame(Strain=est.wet$Strain,IsoFreq=est.wet$PWC,IsoEnv='Wet',
                             Matched_PQ=est.wet$emmean,MisMatched_PQ=est.dry$emmean))

matched_model <- glmer(IsoFreq~Matched_PQ +(1|Strain),mixed.df,family='poisson')
mismatched_model <- glmer(IsoFreq~MisMatched_PQ +(1|Strain),mixed.df,family='poisson')

#qualitative comparison in models
Anova(matched_model)
Anova(mismatched_model)

r.squaredGLMM(matched_model)
r.squaredGLMM(mismatched_model)

#compare models using Akaike Weights (Probability of being best)

aic_vals <- c(AIC(matched_model), AIC(mismatched_model))
delta_aic <- aic_vals - min(aic_vals)
likelihood <- exp(-0.5 * delta_aic)
weights <- likelihood / sum(likelihood)
weights
