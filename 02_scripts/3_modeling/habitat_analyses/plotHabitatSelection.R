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
load(outpath('stepSelectionParamsHSF_DF.rdata'))


# ******************************************************************************
#                                   FUNCTIONS
# ******************************************************************************
plotByFactor <- function(fac, dat) {
  dat$CAT <- dat[,fac]
  ntot <- dat %>%
    group_by(case_, CAT, country) %>%
    dplyr::summarize(tot=n())
  use_vs_avail <- dat %>%
    group_by(case_, CAT, country, newstruc) %>%
    dplyr::summarize(n=n()) %>%
    left_join(ntot, by=c('case_', 'CAT', 'country')) %>%
    mutate(prop = n / tot,
           label = paste0(round(prop * 100, 1), "%"))
  use_vs_avail %>%
    ggplot(aes(x=CAT, y=prop, label=label,
               fill = case_, group = case_)) +
    geom_col(position = position_dodge2()) +
    facet_wrap(~country+newstruc, nrow=2)+
    geom_text(size = 4, vjust = -0.25, position = position_dodge(width = 1)) +
    labs(x = fac, y = "Proportion")+
    scale_fill_brewer(palette = "Paired",
                      breaks=c("FALSE", "TRUE"),
                      labels=c("Available", "Used")) +
    theme_light() +
    guides(fill='none') +
    # coord_flip()+
    big.theme
}


# ******************************************************************************
#                             BRIEF EDA - space use
# ******************************************************************************
data <- ssf %>% 
  filter(!is.na(evi), !is.na(structure)) %>% 
  mutate(tod=factor(tod),
         newstruc = fct_collapse(structure,
                                 other = c('cropland', 'wetland', 'closed bush')))

# plot space use by factor
plotByFactor('season', data)
plotByFactor('tod', data)

# plot available habitat by country
ntot <- data %>%
  group_by(country) %>%
  dplyr::summarize(tot=n())
data %>%
  filter(!case_) %>% 
  group_by(country, structure) %>%
  dplyr::summarize(n=n()) %>%
  left_join(ntot, by=c('country')) %>%
  mutate(prop = n / tot,
         label = paste0(round(prop * 100, 1), "%")) %>% 
  ggplot(aes(x=structure, y=prop, label=label,
             fill = country, group = country)) +
  geom_col(position = position_dodge2()) +
  geom_text(size = 4, vjust = -0.25, position = position_dodge(width = 1)) +
  labs(x = '', y = "Proportion of available habitat")+
  # scale_fill_manual(values=c('#5dd6f5', '#18ab3f')) +
  scale_fill_manual(values=c('gray', 'red')) +
  theme_light() +
  big.theme + theme(text=element_text(size=20))


# plot use/avail side by side
ntot <- data %>%
  group_by(case_, country) %>%
  dplyr::summarize(tot=n())
data %>%
  group_by(case_, country, structure) %>%
  dplyr::summarize(n=n()) %>%
  left_join(ntot, by=c('case_', 'country')) %>%
  mutate(prop = n / tot,
         label = paste0(round(prop * 100, 1), "%")) %>% 
  ggplot(aes(x=structure, y=prop, label=label,
             fill = case_, group = case_)) +
  geom_col(position = position_dodge2()) +
  facet_wrap(~country, nrow=2) +
  geom_text(size = 4, vjust = -0.25, position = position_dodge(width = 1)) +
  labs(x = 'veg structure', y = "Proportion")+
  scale_fill_manual(values=c('darkgray', 'black'),
                    breaks=c("FALSE", "TRUE"),
                    labels=c("Available", "Used")) +
  theme_light() +
  big.theme
