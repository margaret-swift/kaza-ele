# usedHabitatCalibration.R
# Created 03 Sept 2024
# Margaret Swift <margaret.swift@cornell.edu>

# Part of simulation validation portion of project

# paper: https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.03123
# #  Fieberg et al 2018
# 
# 



# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
# packages and HERE
pacman::p_load(here, hmmSSF, momentuHMM)
i_am('02_scripts/5_validation/usedHabitatCalibration.R')

# utilities
source(here('02_scripts', 'utilities.R'))

# loading data
quickload('elephant')