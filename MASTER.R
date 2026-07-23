# MASTER
#
# This MASTER script will conduct the analysis described in the paper 
# "Sensitivity analyses for missing data: A practical guide for planning, 
#  conducting and reporting."
# 
# Written by R Mainzer

# Load packages
library(boot)
library(dplyr)
library(flextable)
library(forcats)
library(gridExtra)
library(ggdist)
library(ggplot2)
library(ggpubr)
library(gtsummary)
library(ggtext)
library(haven)
library(here)
library(mice)
library(plyr)
library(RColorBrewer)
library(scales)
library(StackImpute)
library(stringr)

# Set working directory
setwd(here::here())

# Load functions 
source("functions/Jackknife_Variance.R") # overwrite StackImpute::Jackknife_Variance with RM edited function
source("functions/gcomp_funs.R")
source("functions/gcomp_calcweights.R")
source("functions/cluster_boot_funs.R")

# Clean and describe LSAC data
source("descr_LSAC.R", echo = TRUE)

# Create Figure 2 for paper
source("cr_figure2_v2.R", echo = TRUE)

# Run analysis of continuous outcome
source("sens_params_cont.R", echo = TRUE)
source("an_LSAC_cont.R")

# Run analysis of binary outcome
source("sens_params_bin.R", echo = TRUE)
source("an_LSAC_bin.R")                        

# Run analysis of continuous outcome with interaction
source("an_LSAC_cont_int.R", echo = TRUE)
