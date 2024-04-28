# Non-linearity-and-variablitity-vertebrate-populations

Code and data for the following paper: 'Non-linearity and temporal variability are overlooked components of global vertebrate population dynamics'.

# Data

## Living Planet Database
The population time series we analysed come from the Living Planet Database (public version from 2022 *"LPD2022_public.csv"* downloaded from https://www.livingplanetindex.org/data_portal)

## IUCN Red List 
We used the IUCN Red List status from the 2022 public version of the database (*"IUCNredist.csv"*, downloaded from https://www.iucnredlist.org/)

# Scripts ("code" folder)

+ 00_preliminary_functions.R : 
  **Contains classification functions adapted from Rigal et al. framework (2020).**

+ 01_trajectories_and_variability.R : 
**Calculates population non-linear trajectories (following Rigal et al. framework (2020)) and temporal variability (using the coefficient of variation, the mean squared error, and the consecutive disparity index from Fernandez-Martinez et al. (2018)).**

+ 02_analyses.R : 
**Tests the sources of heterogeneities in populations' non-linearity and temporal variability.**

+ 03_figures.R : 
**Scripts for all figures (including supplementary) in the paper.**

# Figures ("outputs" folder)

Contains all the figures generated from the code.

# References 
Rigal S, Devictor V, Dakos V (2020) A method for classifying and comparing non-linear trajectories of ecological variables. Ecological Indicators 112: 106113

Fernández-Martínez M, Vicca S, Janssens IA, Carnicer J, Martín-Vide J, Peñuelas J (2018) The consecutive disparity index, D: a measure of temporal variability in ecological studies. Ecosphere 9: e02527
