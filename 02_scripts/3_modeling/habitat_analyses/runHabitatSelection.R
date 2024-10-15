# runHabitatSelection.R
# Created 13 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# This file runs through conditional habitat selection model selection to find the
#   most parsimonious resource selection function for elephant habitat use.
# 
# amt vignette: 
#   https://cran.r-project.org/web/packages/amt/vignettes/p4_SSF.html
# A ‘How to’ guide for interpreting parameters in habitat-selection analyses (Fieberg et al 2021): 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2656.13441
# appendix a - hsf examples (from Fieberg et al 2021):
#   In this appendix, we will walk the user through the process of fitting and 
#   interpreting habitat-selection functions.
#   https://conservancy.umn.edu/server/api/core/bitstreams/db16ca0b-ff6a-4710-8a09-2f1aa32d82cd/content
# appendix b - ssf examples (from Fieberg et al 2021):
#   In this appendix, we will walk the user through the process of fitting and 
#   interpreting parameters and output when conducting an iSSF.
#   https://conservancy.umn.edu/server/api/core/bitstreams/63727072-87b1-4b35-b81c-8fd31b8f1e57/content

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/habitat_analyses/runHabitatSelection.R')
source(here::here('02_scripts','utilities.R'))
library(tidyverse)
library(survival) # conditional logistic modeling with clogit()
quickload() #loads common spatial features

# load step selection dat
setOutPath(c('habitat_selection', 'ssf_data'))
# load(outpath('fitmodel.rdata'))
load(outpath('stepSelectionParamsHSF_DF.rdata'))

# landcover metadata
setDataPaths('landcover')
lands.meta <- read.csv(metapath('landcover_meta_2005.csv'))



## REAL QUICK FOR DNRE: find distance from waterholes
setDataPaths('waters')
nat.rast <- terra::rast(rawpath('nat_rast_2017', 'waterdist_evi_raster_dry_2017.tif'))$waterDist
ext = terra::ext(nat.rast)
x = ssf %>% st_as_sf(coords=c('x2_', 'y2_'), crs=32734)
# poly <- data.frame(
#   x=c(ext[1], ext[1], ext[2], ext[2], ext[1]),
#   y=c(ext[3], ext[4], ext[4], ext[3], ext[3])
# ) %>%
#   st_as_sf(coords = c("x", "y"), crs = 32734) %>%
#   summarise(geometry = st_combine(geometry)) %>%
#   st_cast("POLYGON")
# ssf.sub <- st_intersection(x, poly)
water.dists <- terra::extract(nat.rast, ssf[,c(2,4)])
setOutPath(c('habitat_selection', 'ssf_data'))
save(water.dists, file=outpath('waterdists.rdata'))

# ******************************************************************************
#                             RUN CONDITIONAL HSF GLM
# ******************************************************************************

# We are using the clogit() function from the "survival" package, where case 1
# is the step actually taken and case 0 is steps not taken. 

# CONDITIONAL HSF
hsfFunc <- function(f) {
  f <- paste0("case_ ~ ", f, " + strata(inx)")
  message('model == ', f)
  f <- formula(f)
  fit <- clogit(f,
         data=data,
         weights=w,
         method="approximate",
         model=TRUE)
  fit <- cleanUpModel(fit)
  gc()
  fit
}

# regular SSF
m0 <- hsfFunc('country*structure')
gc()
m1 <- hsfFunc('evi')
gc()
m2 <- hsfFunc('evi*season')
gc()
m3 <- hsfFunc('evi*season + structure')
gc()
m4 <- hsfFunc('evi*tod + evi*season + structure')
gc()
m5 <- hsfFunc('evi*season + season*structure')
gc()
mlist.ssf <- list(m1, m2, m3, m4, m5)

# # integrated SSF (using ta and sl)
# m5 <- hsfFunc('evi + structure + sl_ + log_sl + cos_ta')
# gc()
# m6 <- hsfFunc('evi + tod*structure + sl_ + log_sl*structure + cos_ta*structure')
# gc()
# mlist.issf <- list(m5, m6)

# ******************************************************************************
#                            COMPARE MODELS WITH AIC
# ******************************************************************************

# GET AIC AND ARRANGE NICELY
compareModels <- function(mlist) {
  calls <- sapply(mlist, function(e) as.character(e$formula[3]))
  df <- data.frame(model=paste(calls, sep=" : "),
                   AIC=sapply(mlist,  AIC), 
                   LRT=sapply(mlist,  function(e) summary(e)$logtest['test']),
                   LRTp=sapply(mlist, function(e) summary(e)$logtest['pvalue']),
                   rsq= sapply(mlist, function(e) summary(e)$rsq[[1]]))
  df <- arrange(df, AIC)
  df$dAIC <- round(df$AIC - min(df$AIC))
  df
}
aic.df <- compareModels(mlist.ssf)
print(aic.df)

## comparing coefficients for different models
# GETTING MODEL COEFFICIENTS
getCoefs <- function(hsf) {
  tab <- summary(hsf)$coefficients
  vals <- tab[,2]
  se <- tab[,3]
  pvals <- tab[,5]
  star <- rep('', length(pvals))
  star[pvals < 0.1] <- '.'
  star[pvals < 0.05] <- '*'
  star[pvals < 0.01] <- '**'
  star[pvals < 0.001] <- '***'
  df <- data.frame(value=vals, se=se, sig=star)
  df
}

# ******************************************************************************
#                               ANALYZE RESULTS
# ******************************************************************************

# we are going with state x evi and water distance
hsf <- hsfFunc('evi*season + season*structure')
summary(hsf)
coefs <- getCoefs(hsf)
print(coefs)
save(hsf, file=outpath('fitmodelBN.rdata'))


# ******************************************************************************
#                               SAVE RESULTS
# ******************************************************************************


# I think that hsf3 is the best choice, despite having higher AIC than the 
#   model with a state:cover_class interaction. It's less complicated and 
#   honestly the significance values aren't wowing me.
hsf <- hsf3
save(hsf, file=here(outdir, 'habitat_selection', 'hsf.rdata'))


