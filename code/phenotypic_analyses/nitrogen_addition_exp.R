#here we examine the impact of  water history on rhizobium partner quality, and how that interacts with the contemporary nitrogen 
#we use a small subset of strains from the larger experiment. These strains were evolved with a plant and under a low nitrogen environment
#and either came from historic dry or historic wet treatments 
#we examine multiple plant traits of interest, as well as using multivariate approaches (pc for 
#variable compression and rda for full multivariate analysis). Within the Ricks et al manuscript, 
#these analyses correspond to Figure S5

library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(vegan)
library(multcomp)
library(permute)
library(dplyr)


new.d <- read.csv('path/to/plant/phenotypic/nitrogen_addition_data.csv')

#rep was our blocking factor-transform into factor
new.d$Rep <- as.factor(new.d$Rep)
new.d$MergeRep <- paste(new.d$Water,new.d$Rep)


#create root-shoot ratio variable
new.d$RS <- new.d$Belowground_biomass/new.d$Aboveground_biomass


#remove sapmles that had problems-this primarily remove those that had died halfway through the experiment 
#while we measured some data on these, a plant biomass of plant that died 3 weeks before harvest, for example,
#is not an apt comparison
#we note that in this study there were a large amount of plant deaths due to a heat wave while conducting the experiment
new.d <- new.d[new.d$Notes=="",]



#scale variables of interest between 0 and 1
new.d$Aboveground_biomass <- (new.d$Aboveground_biomass-min(new.d$Aboveground_biomass,na.rm=T))/(max(new.d$Aboveground_biomass,na.rm=T)-min(new.d$Aboveground_biomass,na.rm=T))
new.d$Height <- (new.d$Height-min(new.d$Height,na.rm=T))/(max(new.d$Height,na.rm=T)-min(new.d$Height,na.rm=T))
new.d$Leaf_Count <- (new.d$Leaf_Count-min(new.d$Leaf_Count,na.rm=T))/(max(new.d$Leaf_Count,na.rm=T)-min(new.d$Leaf_Count,na.rm=T))
new.d$Belowground_biomass <- (new.d$Belowground_biomass-min(new.d$Belowground_biomass,na.rm=T))/(max(new.d$Belowground_biomass,na.rm=T)-min(new.d$Belowground_biomass,na.rm=T))

#remove missing values, necessary for PCA
new.d <- new.d[!is.na(new.d$Aboveground_biomass) & !is.na(new.d$Height) & 
                 !is.na(new.d$Leaf_Count) &  !is.na(new.d$Belowground_biomass),]



pc.axis <- (prcomp(cbind(new.d$Leaf_Count,new.d$Aboveground_biomass,new.d$Height,new.d$Belowground_biomass),scale=T))
new.d$pc1 <- pc.axis$x[,1]

#variance explained by pc axes
summary(pc.axis)
summary(pc.axis)$importance[2,]

#pc1 is positively correlated with all traits
plot(pc1~Aboveground_biomass,new.d)
plot(pc1~Leaf_Count,new.d)
plot(pc1~Height,new.d)


#run models on pc as well as individual traits
anova(lmer(pc1~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Wet' & new.d$Rhizobia_Strain!='Sterile',]))
anova(lmer(pc1~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[ new.d$Contemporary_water=='Dry' & new.d$Rhizobia_Strain!='Sterile',]))



anova(lmer(Aboveground_biomass~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Wet' & new.d$Rhizobia_Strain!='Sterile',]))
anova(lmer(Belowground_biomass~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Wet' & new.d$Rhizobia_Strain!='Sterile',]))
anova(lmer(Height~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Wet' & new.d$Rhizobia_Strain!='Sterile',]))
anova(lmer(Leaf_Count~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Wet' & new.d$Rhizobia_Strain!='Sterile',]))

anova(lmer(Aboveground_biomass~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Dry' & new.d$Rhizobia_Strain!='Sterile',]))
anova(lmer(Belowground_biomass~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Dry' & new.d$Rhizobia_Strain!='Sterile',]))
anova(lmer(Height~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Dry' & new.d$Rhizobia_Strain!='Sterile',]))
anova(lmer(Leaf_Count~Historic_water*Contemporary_nitrogen +(1|Rep) + (1|Rhizobia_Strain),new.d[new.d$Contemporary_water=='Dry' & new.d$Rhizobia_Strain!='Sterile',]))



#additional multivariate models using rda 
#wet environment
data.sub.wet <- new.d[new.d$Contemporary_water=='Wet' & new.d$Rhizobia_Strain!='Sterile',]
traits_matrix <- data.sub.wet[,colnames(data.sub.wet) %in% c('Leaf_Count','Height','Aboveground_biomass','Belowground_biomass')]
anova(rda(traits_matrix ~ Historic_water*Contemporary_nitrogen + Condition(Rep) ,data=data.sub.wet),by='term')


#dry environment
data.sub.dry <- new.d[new.d$Contemporary_water=='Dry' & new.d$Rhizobia_Strain!='Sterile',]
traits_matrix <- data.sub.dry[,colnames(data.sub.dry) %in% c('Leaf_Count','Height','Aboveground_biomass','Belowground_biomass')]
anova(rda(traits_matrix ~ Historic_water*Contemporary_nitrogen + Condition(Rep) ,data=data.sub.dry),by='term')



