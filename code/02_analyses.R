################################################################################
### Non-linearity and temporal variability are overlooked components of 
### global vertebrate population dynamics
### 04/2024

#### Part 2 - Analyses of population trajectories and variability

################################################################################
# Libraries ---------------------------------------------------------------
library(DHARMa)
library(lme4)
library(dplyr)
library(emmeans)
library(multcomp)
library(reshape2)
library(LMERConvenienceFunctions)
library(glmmTMB)
library(bbmle)
library(effects)
library(ggplot2)

################################################################################
# Load the data -----------------------------------------------------------
save <- read.table("data/non_linear_models.csv",header = T,sep = ",")[,-1]
save <- save %>%
  mutate(LoNL = as.factor(ifelse(acceleration=="constant","linear","non linear")))
d <- save %>% filter(., redlistCategory!="Extinct in the Wild",
                     taxa!="Invertebrates") 
d$Region<-factor(d$Region)
d$System<-as.factor(d$System)
d$taxa<-as.factor(d$taxa)
d$redlistCategory<-factor(d$redlistCategory,
                          levels=c("Critically Endangered",
                                   "Endangered",
                                   "Vulnerable",
                                   "Near Threatened",
                                   "Least Concern",
                                   "Data Deficient"))

################################################################################
############################## NON LINEARITY ###################################
################################################################################

### 1 - Fit the GLMM for non-linearity ------------------------
# GLMM with a binomial error structure and a logit link function
fit_mod <- glmmTMB(LoNL~points+Region+System+taxa+redlistCategory+(1|Binomial),
                   data=d,family=binomial)

### 2 - Test assumptions ---------------------------------------
res_fit<-simulateResiduals(fit_mod)
plot(res_fit) # ok

### 3 - Inference: test anova ----------------------------------
car::Anova(fit_mod) # Points, Region and Taxonomic group effects

### 4 - Estimated Marginal Means -------------------------------
(EMM_hab <- emmeans(fit_mod, ~ System, type="response"))
(EMM_reg <- emmeans(fit_mod, ~ Region, type="response"))
(EMM_taxa <- emmeans(fit_mod, ~ taxa, type="response"))
(EMM_iucn <- emmeans(fit_mod, ~ redlistCategory, type="response"))

### 5 - Pairwise comparisons and CLD --------------------------

# Habitat type : Not really useful as there is no Habitat type effect (see Anova)
# mc_hab <- glht(fit_mod, linfct = mcp(System = "Tukey"))
# tuk.cld <- cld(mc_hab)
# hab.let <- tuk.cld$mcletters$Letters
# hab.let_df <- data.frame(System=levels(d$System),letters=hab.let)

# Region
mc_reg <- glht(fit_mod, linfct = mcp(Region = "Tukey"))
summary(mc_reg)
tuk.cld <- cld(mc_reg) 
reg.let <- tuk.cld$mcletters$Letters
reg.let_df <- data.frame(Region=levels(d$Region),letters=reg.let)

# Taxonomic group
d$taxa <- factor(d$taxa,exclude=NULL)
mc_taxa <- glht(fit_mod, linfct = mcp(taxa = "Tukey"))
summary(mc_taxa)
tuk.cld <- cld(mc_taxa) 
taxa.let <- tuk.cld$mcletters$Letters
taxa.let_df <- data.frame(taxa=levels(d$taxa),letters=taxa.let)

# RLC : Not really useful as there is no Habitat type effect (see Anova)
# mc_rlc <- glht(fit_mod, linfct = mcp(redlistCategory = "Tukey"))
# tuk.cld <- cld(mc_rlc)
# rlc.let <- tuk.cld$mcletters$Letters
# rlc.let_df <- data.frame(redlistCategory=levels(d$redlistCategory),letters=rlc.let)

################################################################################
############################ TEMPORAL VARIABILITY ##############################
################################################################################

### Data adjustments prior to modelling
d$shape_class<-as.factor(d$shape_class)

### 1 - Fit the GLMMs for temporal variability ------------------------
# GLMM with a gamma error structure and a log link function

# VAR 1 : MSE ~ Hab, Reg, Tax, RLC
var_mod <- glmmTMB(MSE~points+Region+System+taxa+redlistCategory+(1|Binomial),
                   data=d,family=Gamma(link="log"))

# VAR 2 : MSE ~ Hab, Reg, Tax, RLC, Trajectory
var_mod_full <- glmmTMB(MSE~points+Region+System+taxa+redlistCategory+shape_class+(1|Binomial),
                    data=d,family=Gamma(link="log"))

# VAR 3 : Using D /!\ NEED TO DO IT ON RAW VALUES !!!!!
var_mod_D <- glmmTMB(D~points+Region+System+taxa+redlistCategory+shape_class+(1|Binomial),
                   data=d,family=Gamma(link="log"))

# VAR 4 : Using CV
var_mod_CV <- glmmTMB(CV~points+Region+System+taxa+redlistCategory+shape_class+(1|Binomial),
                   data=d,family=Gamma(link="log"))

### 2 - Test assumptions ---------------------------------------
plot(simulateResiduals(var_mod)) # VAR1
plot(simulateResiduals(var_mod_full)) # VAR2
plot(simulateResiduals(var_mod_D)) # VAR3
plot(simulateResiduals(var_mod_CV)) # VAR4

### 3 - Inference: test anova ----------------------------------
car::Anova(var_mod) # VAR1
car::Anova(var_mod_full) # VAR2
car::Anova(var_mod_D) # VAR3
car::Anova(var_mod_CV) # VAR4

### 4 - Estimated Marginal Means -------------------------------

# VAR2 - With all trajectory types
EMM_full <- emmeans(var_mod_full, ~ shape_class, type="response")
EMM_full_hab <-  emmeans(var_mod_full, ~ System, type="response")
EMM_full_taxa <-  emmeans(var_mod_full, ~ taxa, type="response")
EMM_full_reg <-  emmeans(var_mod_full, ~ Region, type="response")
EMM_full_iucn <- emmeans(var_mod_full, ~ redlistCategory, type="response")

# VAR3 - With D
EMM_full_D <- emmeans(var_mod_D, ~ shape_class, type="response")
EMM_hab_D <- emmeans(var_mod_D, ~ System, type="response")
EMM_reg_D <- emmeans(var_mod_D, ~ Region, type="response") 
EMM_taxa_D <- emmeans(var_mod_D, ~ taxa, type="response")
EMM_iucn_D <- emmeans(var_mod_D, ~ redlistCategory, type="response")

# VAR4 - With CV
EMM_full_CV <- emmeans(var_mod_CV, ~ shape_class, type="response")
EMM_hab_CV <- emmeans(var_mod_CV, ~ System, type="response")
EMM_reg_CV <- emmeans(var_mod_CV, ~ Region, type="response") 
EMM_taxa_CV <- emmeans(var_mod_CV, ~ taxa, type="response")
EMM_iucn_CV <- emmeans(var_mod_CV, ~ redlistCategory, type="response")

### 5 - Pairwise comparisons and CLD --------------------------

# MAIN TEXT - ONLY MSE (see D and CV below)
# Trajectory type
mc_traj <- glht(var_mod_full, linfct = mcp(shape_class = "Tukey"))
summary(mc_traj)
tuk.cld <- cld(mc_traj) 
traj.let <- tuk.cld$mcletters$Letters
traj.let_df <- data.frame(shape_class=levels(d$shape_class),letters=traj.let)

# Habitat type
mc_hab <- glht(var_mod_full, linfct = mcp(System = "Tukey"))
summary(mc_hab)
tuk.cld <- cld(mc_hab) 
hab.let_full <- tuk.cld$mcletters$Letters
hab.let_full_df <- data.frame(System=levels(d$System),letters=hab.let_full)

# Region
mc_reg <- glht(var_mod_full, linfct = mcp(Region = "Tukey"))
summary(mc_reg)
tuk.cld <- cld(mc_reg) 
reg.let_full <- tuk.cld$mcletters$Letters
reg.let_full_df <- data.frame(Region=levels(d$Region),letters=reg.let_full)

# Taxonomic group
mc_taxa <- glht(var_mod_full, linfct = mcp(taxa = "Tukey"))
summary(mc_taxa)
tuk.cld <- cld(mc_taxa) 
taxa.let_full <- tuk.cld$mcletters$Letters
taxa.let_full_df <- data.frame(taxa=levels(d$taxa),letters=taxa.let_full)


################################################################################
# Same analyses for CV and D 


# D index -----------------------------------------------------------------
# Trajectory type
mc_traj <- glht(var_mod_D, linfct = mcp(shape_class = "Tukey"))
summary(mc_traj)
tuk.cld <- cld(mc_traj) 
traj.let <- tuk.cld$mcletters$Letters
traj.let_D_df <- data.frame(shape_class=levels(d$shape_class),letters=traj.let)

# Habitat type
mc_hab <- glht(var_mod_D, linfct = mcp(System = "Tukey"))
summary(mc_hab)
tuk.cld <- cld(mc_hab) 
hab.let <- tuk.cld$mcletters$Letters
hab.let_D_df <- cbind(data.frame(EMM_hab_D),letters=hab.let)

# Region
mc_reg <- glht(var_mod_D, linfct = mcp(Region = "Tukey"))
summary(mc_reg)
tuk.cld <- cld(mc_reg) 
reg.let <- tuk.cld$mcletters$Letters
reg.let_D_df <- data.frame(data.frame(EMM_reg_D),letters=reg.let)

# Taxonomic group
mc_taxa <- glht(var_mod_D, linfct = mcp(taxa = "Tukey"))
summary(mc_taxa)
tuk.cld <- cld(mc_taxa) 
taxa.let <- tuk.cld$mcletters$Letters
taxa.let_D_df <- data.frame(data.frame(EMM_taxa_D),letters=taxa.let)

# Coefficient of variation ------------------------------------------------
# Trajectory type
mc_traj <- glht(var_mod_CV, linfct = mcp(shape_class = "Tukey"))
summary(mc_traj)
tuk.cld <- cld(mc_traj) 
traj.let <- tuk.cld$mcletters$Letters
traj.let_CV_df <- data.frame(shape_class=levels(d$shape_class),letters=traj.let)

# Habitat type
mc_hab <- glht(var_mod_CV, linfct = mcp(System = "Tukey"))
summary(mc_hab)
tuk.cld <- cld(mc_hab) 
hab.let <- tuk.cld$mcletters$Letters
hab.let_CV_df <- cbind(data.frame(EMM_hab_CV),letters=hab.let)

# Region
mc_reg <- glht(var_mod_CV, linfct = mcp(Region = "Tukey"))
summary(mc_reg)
tuk.cld <- cld(mc_reg) 
reg.let <- tuk.cld$mcletters$Letters
reg.let_CV_df <- data.frame(data.frame(EMM_reg_CV),letters=reg.let)

# Taxonomic group
mc_taxa <- glht(var_mod_CV, linfct = mcp(taxa = "Tukey"))
summary(mc_taxa)
tuk.cld <- cld(mc_taxa) 
taxa.let <- tuk.cld$mcletters$Letters
taxa.let_CV_df <- data.frame(data.frame(EMM_taxa_CV),letters=taxa.let)
