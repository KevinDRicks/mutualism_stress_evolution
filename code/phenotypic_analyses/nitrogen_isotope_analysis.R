#here we examine the impact of  water & nitrogen on rhizobium partner quality,
#specifically in relation to plant nitrogen content. we look at % nitrogen and delta N 15,
#a common N isotope enriched through fixation processes
#Within the Ricks et al manuscript these analyses correspond to Figure S4

library(lme4)
library(lmerTest)
library(emmeans)

data<- read.csv('G:/My Drive/Work/UIUC/Projects/Rhizobia/Drought_evolution/Github/data/greenhouse/nitrogen_isotope_data.csv')

#code ancestor variable as factor
data$ancestor <- as.character(data$ancestor)

#significant variation in nitrogen content dependent on which ancestor these evolved strains map back to
boxplot(N_percent~ancestor,data)
anova(lm(N_percent~ancestor,data))

boxplot(Delta_N~ancestor,data)
anova(lm(Delta_N~ancestor,data))

#weak, and marginal effect of historic nitrogen on plant nitrogen content
n_per_mod <- (lmer(N_percent~Historic_water*Historic_nitrogen + (1|ancestor),data))
anova(n_per_mod)
summary(n_per_mod)
plot(emmeans(n_per_mod,~Historic_nitrogen))
plot(emmeans(n_per_mod,~Historic_nitrogen*Historic_water))

#no impact on delta nitrogen 15 content
n_15_mod <- (lmer(Delta_N~Historic_water*Historic_nitrogen + (1|ancestor),data))
anova(n_15_mod)
summary(n_15_mod)
plot(emmeans(n_15_mod,~Historic_nitrogen*Historic_water))


