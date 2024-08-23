# runHabitatSelection.R
# Created 13 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# This file analyzes a conditional habitat selection model for elephant habitat use.
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
here::i_am('02_scripts/3_modeling/habitat_analyses/analyzeHabitatSelection.R')
source(here::here('02_scripts','utilities.R'))
library(tidyverse)
library(survival) # conditional logistic modeling with clogit()
quickload() #loads common spatial features

# load step selection dat
setOutPath(c('habitat_selection', 'ssf_data'))
load(outpath('fitmodel.rdata'))


# ******************************************************************************
#                                   FUNCTIONS
# ******************************************************************************

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

..glm.ratio <- function(GLM.RESULT, DIGITS = 2, P.DIGITS = 3, CONF.LEVEL = 0.95) {
  ## Extract coefficients and confidence interval
  TABLE     <- cbind(coef = stats::coef(GLM.RESULT), 
                     suppressMessages(stats::confint(GLM.RESULT, level = CONF.LEVEL)))
  ## Turn them into OR
  TABLE.EXP <- round(exp(TABLE), DIGITS)
  colnames(TABLE.EXP)[1] <- "OR"
  ## Extract p-value
  return(cbind(TABLE.EXP, 
               "P" = round(summary(GLM.RESULT)$coef[,4], P.DIGITS))) 
}
# 
# # set up dataframe to fill with predictions based on HSF fit
# predictData <- function(hsf) {
#   # create data
#   evi <- seq(0.05, 0.25, 0.05)
#   structure <- hsf$xlevels$structure
#   season <- hsf$xlevels$season
#   colnames <- c('inx', 'structure', 'season', 'evi', 
#                 'pred', 'pred_lo', 'pred_hi')
#   
#   # set up data frame
#   LE=length(evi); LC=length(structure); LS=length(season)
#   NR=LE*LC*LS; NC=length(colnames);
#   pdata <- data.frame(matrix(0, nrow=NR, ncol=NC))
#   names(pdata) <- colnames
#   
#   # fill pdata
#   pdata$inx <- 1:NR
#   pdata$evi <- evi#rep(evi, by=LS)
#   pdata$season = rep(season, each=LC*LE)
#   pdata$structure = rep(structure, each=LE)
#   
#   # predicted fit
#   fit <- predict(hsf, pdata, type="survival", se.fit=TRUE, interval="confidence")
#   pdata$pred <- exp(fit$fit)
#   pdata$pred_lo <- exp(fit$fit - fit$se.fit)
#   pdata$pred_hi <- exp(fit$fit + fit$se.fit)
#   
#   # return
#   return(pdata)
# }

# ******************************************************************************
#                                 SUMMARIZE HSF FIT
# ******************************************************************************

summary(hsf)
coefs <- getCoefs(hsf)
coefs

# ******************************************************************************
#                                 GRAB ODDS RATIOS
# ******************************************************************************

# get GLM ratios from model
ratios <- ..glm.ratio(hsf) %>% 
  as.data.frame() %>% 
  rownames_to_column(var='covariate') %>% 
  filter(grepl('structure', covariate)) %>% 
  mutate(structure=gsub('.*:|structure','', covariate),
         season=gsub('season|:.*', '', covariate),
         season=ifelse(grepl('structure', season), 'dry_1', season),
         is_base=season=='dry_1') %>% 
  dplyr::select(-covariate)
names(ratios)[1:4] <- c('or', 'lower', 'upper', 'p')
bases <- ratios[ratios$is_base,] %>% 
  dplyr::select(structure,or) %>% 
  rename(base=or)
rat.df = left_join(ratios, bases, by="structure") %>% 
  mutate(or_adj = ifelse(is_base, base, or*base),
         low_adj= ifelse(is_base, lower, lower*base),
         up_adj = ifelse(is_base, upper, upper*base))

# ******************************************************************************
#                                   PLOTTING
# ******************************************************************************


# PLOT ODDS RATIO COMPARISON BY SEASON
rat.df %>% 
  ggplot(aes(x=season, group=season)) + 
  geom_linerange(aes(ymin=low_adj, ymax=up_adj, 
                     color=is_base),
                 linewidth=1,
                 position=position_dodge2(width=0.5)) +
  geom_point(aes(y=or_adj), size=2, 
             position=position_dodge2(width=0.5)) +
  geom_hline(yintercept=1, linetype='dashed', color='darkgray') +
  facet_wrap(~structure, nrow=1) +
  plot.theme + guides(color='none') + 
  scale_color_manual(values=c('gray', 'red')) +
  ylab('odds ratio') + theme(axis.title.x=element_blank())
