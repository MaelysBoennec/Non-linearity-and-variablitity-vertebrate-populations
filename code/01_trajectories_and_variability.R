################################################################################
### Non-linearity and temporal variability are overlooked components of 
### global vertebrate population dynamics
### 04/2024

#### Part 1 - Modeling population trajectories and variability

################################################################################
# Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggpattern)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(tidyr)
library(treemapify)
library(forcats)
library(ggExtra)
library(moments)
library(lme4)
library(caret)
library(multcomp)
library(ggridges)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
source("https://raw.githubusercontent.com/larmarange/JLutils/master/R/lm_right.R")

################################################################################
# Load the data -----------------------------------------------------------

# Get the raw data (Living Planet Database)
lpd <- readr::read_csv(
  here::here("data","LPD2022_public.csv")
) 

# IUCN Red List Categories
RLC <- readr::read_csv(
  here::here("data","IUCNredlist.csv")
) 

# Combine genus and species name to make comparable species names between 
# the two databases
RLC <- RLC %>% 
  unite(col=Binomial, genusName, speciesName, remove=FALSE) %>% 
  dplyr::select(Binomial,redlistCategory,populationTrend)

################################################################################
# Preliminary functions ---------------------------------------------------

# "fun_taxa" transforms the taxonomic Classes originally implemented in the LPD
# into more general taxonomic groups
fun_taxa <- function(L) {
  newL<-list()
  for (i in (1:length(L))){
    if (L[i]=='Aves'){newL<-append(newL,"Birds")}
    if (L[i]=='Mammalia'){newL<-append(newL,"Mammals")}
    if (L[i]=='Reptilia'){newL<-append(newL,"Reptiles")}
    if (L[i]=='Amphibia'){newL<-append(newL,"Amphibians")}
    if (L[i]=='Actinopteri'|L[i]=='Coelacanthi'|L[i]=='Petromyzonti'|L[i]=='Dipneusti'){newL<-append(newL,"Fish")}
    if (L[i]=='Elasmobranchii'|L[i]=='Holocephali'){newL<-append(newL,"Sharks_Rays")}
    if (L[i]=='Myxini'){newL<-append(newL,"Invertebrates")}
  }
  return(as.character(unlist(newL)))
}

# d calculates the consecutive disparity index
# Code from Fernandez-Martinez et al. (2018)
d <- function(ts){
  k <- 0.01 * mean(ts, na.rm=TRUE)
  n <- length(ts)
  ret <- c()
  for(i in 1:(length(ts)-1)){
    ret[i] <- abs(log((ts[i+1] + k) / (ts[i]+k)))
  }
  sum(ret)/(n-1)
}

################################################################################
# Reshape data into long form ---------------------------------------------

lpd.long <- lpd %>%
  tidyr::pivot_longer(cols=c(33:103), names_to = "year", values_to = "abundance") %>% 
  dplyr::mutate(year = as.numeric(year),
                abundance = as.numeric(abundance)) %>% 
  drop_na(abundance) # remove NAs in the abundance column

lpd.long <- lpd.long  %>%
  dplyr::distinct() %>%     # remove duplicate rows
  dplyr::group_by(ID) %>%   # group rows so that each group is one population
  dplyr::mutate(minyear = min(year),
                maxyear = max(year),
                duration = (maxyear - minyear)+1,
                scaleab = (abundance - min(abundance))/(max(abundance) - min(abundance)),
                logab = log(abundance+1),
                taxa = fun_taxa(Class),
                zeros=sum(abundance==0)) %>%
  dplyr::filter(is.finite(scaleab)) %>% 
  ungroup() %>% 
  dplyr::group_by(ID,Binomial,Common_name) %>%
  dplyr::mutate(points=length(year))

# save the data into a csv file
# write.csv(lpd.long,file="data/lpd_long.csv")

################################################################################
# Calculate population change for each population -------------------------

########## LINEAR MODELS ##########
linear.models <- lpd.long %>%
  # Group by the key variables that we want to iterate over
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  # Create a linear model for each group
  # Extract model coefficients using tidy() from the broom package
  dplyr::do(broom::tidy(lm(logab ~ year, .))) %>%
  # Filter out slopes and remove intercept values
  dplyr::filter(term == "year") %>%
  # Get rid of the column term as we don't need it any more
  dplyr::select(-term) %>%
  # Remove groupings
  dplyr::ungroup() 

# Add a column to classify population trends according to the slopes and their
# significance
linear.models <- linear.models %>%
  # sorting the dataframe according the mean trend of taxonomic groups
  dplyr::group_by(taxa) %>%
  dplyr::mutate(mean_trend = mean(estimate)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(taxa = forcats::fct_reorder(taxa, -mean_trend)) %>% 
  dplyr::mutate(classif=
                  dplyr::case_when(
                    estimate>0 & p.value<0.05 ~ "positive",
                    estimate<0 & p.value<0.05 ~ "negative",
                    p.value>=0.05 ~ "no trend"))


########## NON-LINEAR MODELS ##########
# According to Rigal et al. (2020) methodology, the non-linear classification of 
# the time series is possible when time series consist of at least 4 points
nl.sub <- lpd.long[which(lpd.long$points>4),]

# Apply the non linear classification on the loged abundances time series
nl.models <- nl.sub %>% 
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  dplyr::do(class.trajectory(.$logab,.$year)) %>%
  dplyr::ungroup() 


# Calculate the magnitude of the curvature (using the radius of curvature)
# Attribute the usual linear classification (lin_class)
# Attribute "good news" or "bad news" depending on the non-linear trajectories
#     considering declines and concave trajectories are "bad news" and increases and convex trajectories are "good news"
# Re-classify the 9 types of trajectories into 6 more simple classes (dir2)
nl.models <- nl.models %>% 
  dplyr::mutate(magnitude=
                  dplyr::case_when(
                    direction=="increase" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p3+first_order_coefficient)^2)^(3/2))),
                    direction=="decrease" & second_order_pvalue<=0.05 ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p2+first_order_coefficient)^2)^(3/2))),
                    acceleration=="concave"|acceleration=="convex" ~ abs((2*second_order_coefficient)/((1+(2*second_order_coefficient*p1+first_order_coefficient)^2)^(3/2))),
                    second_order_pvalue>0.05 ~ 0),
                lin_class=
                  dplyr::case_when(
                    linear_slope>0 & linear_slope_pvalue<0.05 ~ "positive",
                    linear_slope<0 & linear_slope_pvalue<0.05 ~ "negative",
                    linear_slope_pvalue>=0.05 ~ "no trend"),
                gnbn=dplyr::case_when((direction=="decrease"|acceleration=="concave") ~ "bad",
                                      (direction=="increase"|acceleration=="convex") ~ "good",
                                      shape_class=="stable_constant" ~ "no"),
                dir2=dplyr::case_when((direction=="decrease" & acceleration=="constant")~"decrease linear",
                                      (direction=="decrease" & acceleration!="constant")~"decrease non linear",
                                      (direction=="increase" & acceleration=="constant")~"increase linear",
                                      (direction=="increase" & acceleration!="constant")~"increase non linear",
                                      (direction=="stable" & acceleration=="constant")~"no trend linear",
                                      (direction=="stable" & acceleration!="constant")~"no trend non linear")
  )

# Same for raw abundance
nl.models.raw <- nl.sub %>% 
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  dplyr::do(class.trajectory(.$abundance,.$year)) %>%
  dplyr::ungroup() 


########## SUB-SELECTION MODELS ##########

# We selected population time series (TS) with 20 years of data at least
nl.sub20 <- lpd.long[which(lpd.long$points>=20),]

# linear models # 32,211 TS
length(unique(linear.models$Binomial)) # 4883 unique species

# non linear models # 21,841 TS
# populations w/ 10 years data
final.nl10 <- nl.models %>% filter(points>=10) # 14,695 TS

# populations w/ 20 years data - log abundance
final.nl <- nl.models %>% filter(points>=20)   # 6,437 TS

# populations w/ 20 years data - raw abundance
final.nl.raw <- nl.models.raw %>% filter(points>=20)
final.nl$D<-final.nl.raw$D # otherwise D is double log tranform in final.nl

# Some transformation to make the dataset manageable
final.nl$dir2<-as.factor(final.nl$dir2)
final.nl$System<-as.factor(final.nl$System)
final.nl$taxa<-as.factor(final.nl$taxa)
final.nl$Latitude<-as.numeric(final.nl$Latitude)
final.nl <- final.nl %>% 
  mutate(full_realm=case_when(T_realm!="NULL"~T_realm,
                              FW_realm!="NULL"~FW_realm,
                              M_realm!="NULL"~M_realm),
         LoNL = ifelse(acceleration=="constant","linear","non linear")) 

final.nl<-left_join(final.nl,RLC,by="Binomial")
levels(final.nl$dir2)<-c("Decrease\nlinear","Decrease\nnon linear","Increase\nlinear",
                         "Increase\nnon linear","No trend\nlinear","No trend\nnon linear")

final.nl$redlistCategory<-as.factor(final.nl$redlistCategory) 

# save the data into a csv file
write.csv(final.nl,file="data/non_linear_models.csv")

################################################################################
### First small exploration
length(unique(final.nl$Binomial))  # 1,257 unique species in the final dataset
table(final.nl$shape_class) # Number of TS within each trajectory type
table(final.nl$acceleration=="constant")    # 2,887 non-linear trajectories : 44.8 % (log abundance)
table(final.nl.raw$acceleration=="constant")# 45.8 % (raw abundance)
table(final.nl10$acceleration=="constant") # With 10 years data : 30 %

################################################################################
################################################################################

################################################################################
#### SUPPLEMENTARY ANALYSIS - SM2 - Analysis of the populations trajectories 
#### classification sensitivity to the log-transformation
################################################################################

################################################################################
# Apply the non linear classification with a GLM on raw abundances time series
nl.models.glm <- nl.sub %>% 
  dplyr::group_by(ID,Binomial,Common_name,
                  Location,Country,Region,Latitude,Longitude,
                  Class,Order,System,
                  T_realm,T_biome,FW_realm,FW_biome,M_realm,M_biome,
                  minyear,maxyear,duration,points,
                  taxa) %>%
  dplyr::do(class.trajectory.glm(.$abundance,.$year)) %>%
  dplyr::ungroup() 


# Filter the shape class and ID of each population for each statistical framework
nl.mods.log <- nl.models %>% dplyr::filter(points>=20) %>% 
  dplyr::select(ID,shape_class) %>% 
  dplyr::rename(shape_class_log=shape_class)
# table(nl.mods.log$shape_class_log)

nl.mods.glm <- nl.models.glm %>% dplyr::filter(points>=20) %>% 
  dplyr::select(ID,shape_class) %>% 
  dplyr::rename(shape_class_glm=shape_class)
# table(nl.mods.glm$shape_class_glm)

# Create a data frame to combine both classifications
nl.mods.comp <- left_join(nl.mods.log,nl.mods.glm,by="ID")
nl.mods.comp$shape_class_log<-as.factor(nl.mods.comp$shape_class_log)
nl.mods.comp$shape_class_glm<-as.factor(nl.mods.comp$shape_class_glm)

# Create and plot the confusion matrix
cm <- confusionMatrix(data=nl.mods.comp$shape_class_log, 
                      reference = nl.mods.comp$shape_class_glm)

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

(comp_glm_log <- ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) + theme_bw() +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "LM on log transformed abundances",y = "GLM on raw abundances") +
  theme(axis.text.x = element_text(angle=45,hjust=1),text=element_text(size=15)))

################################################################################
# Do the same analysis for rare population time series

rare.pop.long <- lpd.long %>% filter(zeros>0 & points >=20)
ID_rare_pop <- unique(rare.pop.long$ID)

nl.mods.log.rare <- nl.models %>% dplyr::filter(ID %in% ID_rare_pop) %>% 
  dplyr::select(ID,shape_class) %>% dplyr::rename(shape_class_log=shape_class)

nl.mods.glm.rare <- nl.models.glm %>% dplyr::filter(ID %in% ID_rare_pop) %>% 
  dplyr::select(ID,shape_class) %>% dplyr::rename(shape_class_glm=shape_class)

nl.mods.comp.rare <- left_join(nl.mods.log.rare,nl.mods.glm.rare,by="ID")
nl.mods.comp.rare$shape_class_log<-as.factor(nl.mods.comp.rare$shape_class_log)
nl.mods.comp.rare$shape_class_glm<-as.factor(nl.mods.comp.rare$shape_class_glm)
cm <- confusionMatrix(data=nl.mods.comp.rare$shape_class_log, 
                      reference = nl.mods.comp.rare$shape_class_glm)

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

(comp_glm_log_rare <- ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) + theme_bw() +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "LM on log transformed abundances",y = "GLM on raw abundances") +
  theme(axis.text.x = element_text(angle=45,hjust=1),text=element_text(size=15)))


################################################################################
#### SUPPLEMENTARY ANALYSIS - SM3 - Impact of the number of years sampled and 
#### starting year of the time series on the detection of non-linearity
################################################################################

# Correspondence between the number of points within the analyzed 
# time series and the duration of the time series
ggplot(data=final.nl, aes(x=points,y=duration))+geom_point(alpha=0.2)+
  theme_bw()+xlab("Number of points")+ylab("Duration (number of years)")+
  theme(text=element_text(size=15))

# Sorting into 5 years windows
final.nl$points <- as.numeric(final.nl$points)
final.nl$points5cl <- cut(final.nl$points, seq(15,70,by=5))
final.nl$duration5cl <- cut(final.nl$duration, seq(15,70,by=5))

final.nl <- final.nl %>% 
  mutate(LoNL = as.factor(ifelse(acceleration=="constant","linear","non linear")))

# Testing the effects

# GLM with a binamial error structure
test_dur <- glm(LoNL ~ points + first_X_value, data=final.nl, family="binomial")
plot(simulateResiduals(test_dur)) # assumptions
summary(test_dur) #
car::Anova(test_dur)
