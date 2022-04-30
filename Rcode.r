#Temperature variation magnifies chlorpyrifos toxicity differently between larval and adult mosquitoes
#Vienna Delnat, Tam Tran, Julie Verheyen, Khuong Van Dinh, Lizanne Janssens and Robby Stoks
#Science of the Total Environment (2019)
#R code tested on 17/03/2022

#####Packages#####

install.packages("afex")
install.packages("lme4")
install.packages("car")     
install.packages("lsmeans")
install.packages("effects")     

library(afex)
library(lme4)
library(car)     
library(lsmeans)
library(effects)     

sessionInfo()

##Set working directory to source file
#RStudio -> Session -> Set Working Directory...-> To Source File Location

#####Datasets#####

#Dataset of mortality (binomial)
dataMortality=read.csv("Delnat-et-al_DTV-mortality-binomial.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("DTV", "Pesticide", "Stage", "Group", "Sex", "Replicate", "Day", "Number")
dataMortality[Factors] <- do.call(cbind.data.frame, lapply(dataMortality[Factors], as.factor))
##Set levels in factor
dataMortality$Pesticide=factor(dataMortality$Pesticide,levels=c("Control","Chlorpyrifos")) 
dataMortality$Group=factor(dataMortality$Group, levels=c("larvae", "male", "female"))
str(dataMortality) 
#Subset
dataMortality2d=subset(dataMortality, Mortality2d!="NA")
str(dataMortality2d) 
dataLARVAE=subset(dataMortality2d, Group=="larvae")
str(dataLARVAE) 
dataMALE=subset(dataMortality2d, Group=="male")
str(dataMALE) 
dataFEMALE=subset(dataMortality2d, Group=="female")
str(dataFEMALE) 

#Dataset of CTmax
dataCTmax=read.csv("Delnat-et-al_DTV-CTmax.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("DTV", "Pesticide", "Stage", "Group", "Sex", "Replicate", "Date", "Time")
dataCTmax[Factors] <- do.call(cbind.data.frame, lapply(dataCTmax[Factors], as.factor))
##Set levels in factor
dataCTmax$Pesticide=factor(dataCTmax$Pesticide,levels=c("Control","Chlorpyrifos")) 
dataCTmax$Group=factor(dataCTmax$Group, levels=c("larvae", "male", "female"))
#Subset
dataCTmax=subset(dataCTmax, CTmax !="NA" & MassPerLarvae != "NA")
str(dataCTmax)
dataCTmaxLARVAE=subset(dataCTmax, Group=="larvae")
str(dataCTmaxLARVAE)
dataCTmaxMALE=subset(dataCTmax, Group=="male")
str(dataCTmaxMALE)
dataCTmaxFEMALE=subset(dataCTmax, Group=="female")
str(dataCTmaxFEMALE)
dataCTmaxSolvent=subset(dataCTmax, Pesticide=="Control")
str(dataCTmaxSolvent)


#####Mortality4d - binomial - Pre-pesticide exposure period#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Generalized linear mixed models with a binomial error structure and the logit link
#Correct for pseudoreplication (Replicate = Vial)
glmerMortality4d=glmer(Mortality4d ~ DTV*Group + (1|Replicate), data=dataMortality, na.action=na.omit, family=binomial(link=logit),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerMortality4d, type="III") 

#Quick effects plot - not used in manuscript
plot(effect(mod=glmerMortality4d, term="DTV*Group"))
#Emmeans and standard errors for figure
PreExposureMortalityPlotData <- summary(lsmeans(glmerMortality4d, ~ DTV*Group, type = "response"))

#Assumption - Dispersion parameter
glmMortality4d=glm(Mortality4d ~ DTV*Group, data=dataMortality, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmMortality4d) 


#####Mortality2d - binomial - Pesticide exposure period#####

# #Interaction in model --> use set_sum_contrasts() and type=3 in Anova
# set_sum_contrasts()
# 
# #Generalized linear mixed models with a binomial error structure and the logit link
# #Correct for pseudoreplication (Replicate = Vial)
# glmerMortality2d=glmer(Mortality2d ~ DTV*Pesticide*Group + (1|Replicate), data=dataMortality2d, na.action=na.omit, family=binomial(link=logit),
#                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
# Anova(glmerMortality2d, type="III")
# 
# #Posthoc test - contrast analysis with fdr correction
# pairs(lsmeans(glmerMortality2d, ~Pesticide|Group, adjust="fdr"))
# pairs(lsmeans(glmerMortality2d, ~Group|Pesticide, adjust="fdr"))
# 
# #Quick effects plot - not used in manuscript
# plot(effect(mod=glmerMortality2d, term="DTV*Pesticide*Group"),type = "response")
# #Emmeans and standard errors for figure
# PesticideExposureMortalityPlotData <- summary(lsmeans(glmerMortality2d, ~ DTV*Pesticide*Group, type = "response"))
# 
# #Assumption - Dispersion parameter
# glmMortality2d=glm(Mortality2d ~ DTV*Pesticide*Group, data=dataMortality2d, na.action=na.omit, family=quasibinomial(link=logit))
# summary(glmMortality2d)


#####Mortality2d - binomial - Pesticide exposure period - Larvae#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Generalized linear mixed models with a binomial error structure and the logit link
#Correct for pseudoreplication (Replicate = Vial)
glmerMortalityLARVAE=glmer(Mortality2d ~ DTV*Pesticide + (1|Replicate), data=dataLARVAE, na.action=na.omit, family=binomial(link=logit),
                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerMortalityLARVAE, type="III") 

#Quick effects plot - not used in manuscript
plot(effect(mod=glmerMortalityLARVAE, term="DTV*Pesticide"),type = "response")
#Emmeans and standard errors for figure
MortalityLarvaePlotData <- summary(lsmeans(glmerMortalityLARVAE, ~ DTV*Pesticide, type = "response"))

#Assumption - Dispersion parameter
glmMortalityLARVAE=glm(Mortality2d ~ DTV*Pesticide, data=dataLARVAE, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmMortalityLARVAE) 


#####Mortality2d - binomial - Pesticide exposure period - Male#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Generalized linear mixed models with a binomial error structure and the logit link
#Correct for pseudoreplication (Replicate = Vial)
glmerMortalityMALE=glmer(Mortality2d ~ DTV*Pesticide + (1|Replicate), data=dataMALE, na.action=na.omit, family=binomial(link=logit),
                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerMortalityMALE, type="III") 

#Quick effects plot - not used in manuscript
plot(effect(mod=glmerMortalityMALE, term="DTV*Pesticide"),type = "response")
#Emmeans and standard errors for figure
MortalityMalePlotData <- summary(lsmeans(glmerMortalityMALE, ~ DTV*Pesticide, type = "response"))

#Assumption - Dispersion parameter
glmMortalityMALE=glm(Mortality2d ~ DTV*Pesticide, data=dataMALE, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmMortalityMALE) 


#####Mortality2d - binomial - Pesticide exposure period - Female#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Generalized linear mixed models with a binomial error structure and the logit link
#Correct for pseudoreplication (Replicate = Vial)
glmerMortalityFEMALE=glmer(Mortality2d ~ DTV*Pesticide + (1|Replicate), data=dataFEMALE, na.action=na.omit, family=binomial(link=logit),
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerMortalityFEMALE, type="III") 

#Quick effects plot - not used in manuscript
plot(effect(mod=glmerMortalityFEMALE, term="DTV*Pesticide"),type = "response")
#Emmeans and standard errors for figure
MortalityFemalePlotData <- summary(lsmeans(glmerMortalityFEMALE, ~ DTV*Pesticide, type = "response"))

#Assumption - Dispersion parameter
glmMortalityFEMALE=glm(Mortality2d ~ DTV*Pesticide, data=dataFEMALE, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmMortalityFEMALE) 


#####CTmax#####

# #General linear mixed model with a normal error structure and the identity link
# #Correct for pseudoreplication (Replicate = Vial)
# lmerCTmax = lmer(CTmax ~ DTV*Pesticide*Group + MassPerLarvae + (1|Replicate) + (1|Date), data=dataCTmax) 
# Anova(lmerCTmax, type=3, white.adjust=TRUE)
# 
# #Posthoc test - contrast analysis with fdr correction
# pairs(lsmeans(lmerCTmax, ~Group|DTV*Pesticide, adjust="fdr")) 
# pairs(lsmeans(lmerCTmax, ~DTV*Pesticide|Group, adjust="fdr")) 
# 
# #Quick effects plot - not used in manuscript
# plot(effect(mod=lmerCTmax, term="DTV*Pesticide*Group"))
# #Emmeans and standard errors for figure
# CTmaxLarvaePlotData <- summary(emmeans(lmerCTmax, ~ DTV*Pesticide*Group, type = "response"))
# 
# #Assumption - Normality of residuals
# shapiro.test(resid(lmerCTmax))                  
# hist(resid(lmerCTmax))    
# #Assumption - Homogeneity of variance
# leveneTest(CTmax ~ DTV*Pesticide*Group, data = dataCTmax)
# #Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
# aggregate(CTmax ~ DTV*Pesticide*Group, data = dataCTmax, var)
# 
# #Outliers and influential observations
# outlierTest(lmerCTmax)
# cd=cooks.distance(lmerCTmax); which(cd>1)
# influenceIndexPlot(lmerCTmax, vars = c("studentized", "Bonf"))


#####CTmax - Solvent Control#####

#General linear mixed model with a normal error structure and the identity link
#Correct for pseudoreplication (Replicate = Vial)
lmerCTmaxSolvent = lmer(CTmax ~ DTV*Group + MassPerLarvae + (1|Replicate) + (1|Date), data=dataCTmaxSolvent) 
Anova(lmerCTmaxSolvent, type=3, white.adjust=TRUE)

#Posthoc test - contrast analysis with fdr correction
pairs(lsmeans(lmerCTmaxSolvent, ~Group, adjust="fdr")) 

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerCTmaxSolvent, term="DTV*Group"))
#Emmeans and standard errors for figure
CTmaxSolventPlotData <- summary(emmeans(lmerCTmaxSolvent, ~ DTV*Group, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerCTmaxSolvent))                  
hist(resid(lmerCTmaxSolvent))    
#Assumption - Homogeneity of variance
leveneTest(CTmax ~ DTV*Group, data = dataCTmaxSolvent)
#Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
aggregate(CTmax ~ DTV*Group, data = dataCTmaxSolvent, var)

#Outliers and influential observations
outlierTest(lmerCTmaxSolvent)
cd=cooks.distance(lmerCTmaxSolvent); which(cd>1)
influenceIndexPlot(lmerCTmaxSolvent, vars = c("studentized", "Bonf"))


#####CTmax - Larvae#####

#General linear mixed model with a normal error structure and the identity link
#Correct for pseudoreplication (Replicate = Vial)
lmerCTmaxLarvae = lmer(CTmax ~ DTV*Pesticide + MassPerLarvae + (1|Replicate) + (1|Date), data=dataCTmaxLARVAE) 
Anova(lmerCTmaxLarvae, type=3, white.adjust=TRUE)

#Posthoc test - contrast analysis with fdr correction
pairs(lsmeans(lmerCTmaxLarvae, ~DTV|Pesticide, adjust="fdr")) 
pairs(lsmeans(lmerCTmaxLarvae, ~Pesticide|DTV, adjust="fdr")) 

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerCTmaxLarvae, term="DTV*Pesticide"))
#Emmeans and standard errors for figure
CTmaxLarvaePlotData <- summary(emmeans(lmerCTmaxLarvae, ~ DTV*Pesticide, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerCTmaxLarvae))                  
hist(resid(lmerCTmaxLarvae))    
#Assumption - Homogeneity of variance
leveneTest(CTmax ~ DTV*Pesticide, data = dataCTmaxLARVAE)
#Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
aggregate(CTmax ~ DTV*Pesticide, data = dataCTmaxLARVAE, var)

#Outliers and influential observations
outlierTest(lmerCTmaxLarvae)
cd=cooks.distance(lmerCTmaxLarvae); which(cd>1)
influenceIndexPlot(lmerCTmaxLarvae, vars = c("studentized", "Bonf"))

# #Start the permutation on lmerCTmaxLarvae
# #Method based on: http://www.uvm.edu/~dhowell/StatPages/Permutation%20Anova/PermTestsAnova.html
# ANOVA <- Anova(lmerCTmaxLarvae, type="III") #F and df from ANOVA
# AnovaChi <- getElement(ANOVA, "Chisq")      #Extract F or chi-square values
# Ftemp <-  AnovaChi[2]                       #Save the F value per factor
# Fpest <-  AnovaChi[3]
# Finteract <-  AnovaChi[5]
# 
# # Now start resampling
# nreps <- 5000           #Number of permutations
# FT <- numeric(nreps)    #Set up space to store F values as calculated.
# FP <- numeric(nreps)    #FT = temperature, FP = pesticide, FI = interaction
# FI <- numeric(nreps)
# FT[1] <- Ftemp          #Save first F of our 5000 permutations (model is already run above)
# FP[1] <- Fpest
# FI[1] <- Finteract
# for (i in 2:nreps) {
#   newants <- sample(dataCTmaxLARVAE$CTmax, 447) #number of observations of your data/subset
#   mod2 <- lmer(newants ~ DTV*Pesticide + MassPerLarvae + (1|Replicate) + (1|Date), data=dataCTmaxLARVAE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
#   ANOVA <- Anova(mod2, type="III")
#   AnovaChi <- getElement(ANOVA, "Chisq")
#   FT[i] <- AnovaChi[2]
#   FP[i] <- AnovaChi[3]
#   FI[i] <- AnovaChi[5]
# }
# 
# #The addition of "+ .Machine$double.eps" is an aid against two numbers 
# #that differ only by floating point computer calculations at the extreme.
# probT <- length(FT[FT >= Ftemp + .Machine$double.eps ^0.5])/nreps
# probP <- length(FP[FP >= Fpest + .Machine$double.eps ^0.5])/nreps       
# probI <- length(FI[FI >= Finteract + .Machine$double.eps ^0.5])/nreps
# 
# #p-values of your model (F and df values already run above in ANOVA)
# cat("The probability value for DTV is ",probT, "\n")
# cat("The probability value for pesticide is ", probP, "\n")
# cat("The probability value for DTV-pesticide interaction is ", probI, "\n")


#####CTmax - Male#####

#General linear mixed model with a normal error structure and the identity link
#Correct for pseudoreplication (Replicate = Vial)
lmerCTmaxMale<-lmer(CTmax ~ DTV*Pesticide + MassPerLarvae + (1|Replicate) + (1|Date), data = dataCTmaxMALE)
Anova(lmerCTmaxMale, type=3, White.adjust=TRUE)

#Posthoc test - contrast analysis with fdr correction
pairs(lsmeans(lmerCTmaxLarvae, ~DTV, adjust="fdr")) 

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerCTmaxMale, term="DTV*Pesticide"))
#Emmeans and standard errors for figure
CTmaxMalePlotData <- summary(emmeans(lmerCTmaxMale, ~ DTV*Pesticide, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerCTmaxMale))                  
hist(resid(lmerCTmaxMale))    
#Assumption - Homogeneity of variance
leveneTest(CTmax ~ DTV*Pesticide, data = dataCTmaxMALE)
#Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
aggregate(CTmax ~ DTV*Pesticide, data = dataCTmaxMALE, var)

#Outliers and influential observations
outlierTest(lmerCTmaxMale)
cd=cooks.distance(lmerCTmaxMale); which(cd>1)
influenceIndexPlot(lmerCTmaxMale, vars = c("studentized", "Bonf"))

# #Start the permutation on lmerCTmaxMale
# #Method based on: http://www.uvm.edu/~dhowell/StatPages/Permutation%20Anova/PermTestsAnova.html
# ANOVA <- Anova(lmerCTmaxMale, type="III") #F and df from ANOVA
# AnovaChi <- getElement(ANOVA, "Chisq")    #Extract F or chi-square values
# Ftemp <-  AnovaChi[2]                     #Save the F value per factor
# Fpest <-  AnovaChi[3]
# Finteract <-  AnovaChi[5]
# 
# # Now start resampling
# nreps <- 5000           #Number of permutations
# FT <- numeric(nreps)    #Set up space to store F values as calculated.
# FP <- numeric(nreps)    #FT = temperature, FP = pesticide, FI = interaction
# FI <- numeric(nreps)
# FT[1] <- Ftemp          #Save first F of our 5000 permutations (model is already run above)
# FP[1] <- Fpest
# FI[1] <- Finteract
# for (i in 2:nreps) {
#   newsamples <- sample(dataCTmaxMALE$CTmax, 425) #Number of observations of your data/subset
#   M2_males <- lmer(newsamples ~ DTV*Pesticide + MassPerLarvae + (1|Replicate) + (1|Date), data=dataCTmaxMALE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
#   ANOVA <- Anova(M2_males, type="III")
#   AnovaChi <- getElement(ANOVA, "Chisq")
#   FT[i] <- AnovaChi[2]
#   FP[i] <- AnovaChi[3]
#   FI[i] <- AnovaChi[5]
# }
# 
# #The addition of "+ .Machine$double.eps" is an aid against two numbers 
# #that differ only by floating point computer calculations at the extreme.
# probT <- length(FT[FT >= Ftemp + .Machine$double.eps ^0.5])/nreps
# probP <- length(FP[FP >= Fpest + .Machine$double.eps ^0.5])/nreps       
# probI <- length(FI[FI >= Finteract + .Machine$double.eps ^0.5])/nreps
# 
# #p-values of your model (F and df values already run above in ANOVA)
# cat("The probability value for DTV is ",probT, "\n")
# cat("The probability value for pesticide is ", probP, "\n")
# cat("The probability value for DTV-pesticide interaction is ", probI, "\n")


#####CTmax - Female#####

lmerCTmaxFemale <-lmer(CTmax ~ DTV*Pesticide + MassPerLarvae + (1|Replicate) + (1|Date), data = dataCTmaxFEMALE)
Anova(lmerCTmaxFemale, type=3, White.adjust=TRUE)

#Posthoc test - contrast analysis with fdr correction
pairs(lsmeans(lmerCTmaxFemale, ~DTV|Pesticide, adjust="fdr")) 
pairs(lsmeans(lmerCTmaxFemale, ~Pesticide|DTV, adjust="fdr")) 

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerCTmaxFemale, term="DTV*Pesticide"))
#Emmeans and standard errors for figure
CTmaxFemalePlotData <- summary(emmeans(lmerCTmaxFemale, ~ DTV*Pesticide, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerCTmaxFemale))                  
hist(resid(lmerCTmaxFemale))    
#Assumption - Homogeneity of variance
leveneTest(CTmax ~ DTV*Pesticide, data = dataCTmaxFEMALE)
#Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
aggregate(CTmax ~ DTV*Pesticide, data = dataCTmaxFEMALE, var)

#Outliers and influential observations
outlierTest(lmerCTmaxFemale)
cd=cooks.distance(lmerCTmaxFemale); which(cd>1)
influenceIndexPlot(lmerCTmaxFemale, vars = c("studentized", "Bonf"))

# #Start the permutation on lmerCTmaxFemale
# #Method based on: http://www.uvm.edu/~dhowell/StatPages/Permutation%20Anova/PermTestsAnova.html
# ANOVA_fe <- Anova(lmerCTmaxFemale, type="III") #F and df from ANOVA
# AnovaChi_fe <- getElement(ANOVA_fe, "Chisq")   #Extract F or chi-square values
# Ftemp_fe <-  AnovaChi_fe[2]                    #Save the F value per factor
# Fpest_fe <-  AnovaChi_fe[3]
# Finteract_fe <-  AnovaChi_fe[5]
# 
# # Now start resampling on Females' data
# nreps <- 5000              #Number of permutations
# FT_fe <- numeric(nreps)    #Set up space to store F values as calculated.
# FP_fe <- numeric(nreps)    #FT = temperature, FP = pesticide, FI = interaction
# FI_fe <- numeric(nreps)
# FT_fe[1] <- Ftemp_fe       #Save first F of our 5000 permutations (model is already run above)
# FP_fe[1] <- Fpest_fe
# FI_fe[1] <- Finteract_fe
# for (i in 2:nreps) {
#   newsamples <- sample(dataCTmaxFEMALE$CTmax, 397) #Number of observations of your data/subset
#   M4_females <- lmer(newsamples ~ DTV*Pesticide + MassPerLarvae + (1|Replicate)+ (1|Date),
#                      data=dataCTmaxFEMALE,control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
#   ANOVA_fe <- Anova(M4_females, type="III")
#   AnovaChi_fe <- getElement(ANOVA_fe, "Chisq")
#   FT_fe[i] <- AnovaChi_fe[2]
#   FP_fe[i] <- AnovaChi_fe[3]
#   FI_fe[i] <- AnovaChi_fe[5]
# }
# 
# #The addition of "+ .Machine$double.eps" is an aid against two numbers 
# #that differ only by floating point computer calculations at the extreme.
# probT <- length(FT_fe[FT_fe >= Ftemp_fe + .Machine$double.eps ^0.5])/nreps
# probP <- length(FP_fe[FP_fe >= Fpest_fe + .Machine$double.eps ^0.5])/nreps       
# probI <- length(FI_fe[FI_fe >= Finteract_fe + .Machine$double.eps ^0.5])/nreps
# 
# #p-values of your model (F and df values already run above in ANOVA)
# cat("The probability value for DTV is ",probT, "\n")
# cat("The probability value for pesticide is ", probP, "\n")
# cat("The probability value for DTV-pesticide interaction is ", probI, "\n")       


######Save Rdata######
save.image(file="Chapter4_Rdata_20220318_NotPublished.Rdata")
