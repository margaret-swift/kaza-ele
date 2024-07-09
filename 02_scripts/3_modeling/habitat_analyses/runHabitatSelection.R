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

# ******************************************************************************
#                             BRIEF EDA - space use
# ******************************************************************************
data <- ssf %>% 
  filter(!is.na(evi), !is.na(structure)) %>% 
  mutate(tod=factor(TOD),
         newstruc = fct_collapse(structure,
             other = c('cropland', 'wetland', 'closed bush')))
plotByFactor <- function(fac, dat) {
  dat$CATEGORY <- dat[,fac]
  ntot <- dat %>%
    group_by(case_, CATEGORY) %>%
    dplyr::summarize(tot=n())
  use_vs_avail <- dat %>%
    group_by(case_, CATEGORY, newstruc) %>%
    dplyr::summarize(n=n()) %>%
    left_join(ntot, by=c('case_', 'CATEGORY')) %>%
    mutate(prop = n / tot,
           label = paste0(round(prop * 100, 1), "%"))

  use_vs_avail %>%
    ggplot(aes(x=CATEGORY, y=prop, label=label,
               fill = case_, group = case_)) +
    geom_col(position = position_dodge2()) +
    facet_wrap(~newstruc, nrow=2)+
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

plotByFactor('season', data)
plotByFactor('tod', data)

ntot <- data %>%
  group_by(case_) %>%
  dplyr::summarize(tot=n())
data %>%
  group_by(case_, structure) %>%
  dplyr::summarize(n=n()) %>%
  left_join(ntot, by=c('case_')) %>%
  mutate(prop = n / tot,
         label = paste0(round(prop * 100, 1), "%")) %>% 
  ggplot(aes(x=structure, y=prop, label=label,
             fill = case_, group = case_)) +
  geom_col(position = position_dodge2()) +
  geom_text(size = 4, vjust = -0.25, position = position_dodge(width = 1)) +
  labs(x = 'veg structure', y = "Proportion")+
  scale_fill_manual(values=c('darkgray', 'black'),
                    breaks=c("FALSE", "TRUE"),
                    labels=c("Available", "Used")) +
  theme_light() +
  big.theme

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

names(data) <- tolower(names(data))
# regular SSF
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
hsf <- m5#hsfFunc("evi*tod + structure*cos_ta + structure*log_sl")
summary(hsf)
coefs <- getCoefs(hsf)
print(coefs)
save(hsf, file=outpath('fitmodel.rdata'))

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
         low_adj = ifelse(is_base, lower, lower*base),
         up_adj = ifelse(is_base, upper, upper*base))

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

# set up dataframe to fill with predictions
predictData <- function(id, hsf) {
  evi <- seq(0.05, 0.25, 0.05)
  water_dist <- mean(data$water_dist, na.rm=TRUE)
  water_dist <- mean(data$water_dist, na.rm=TRUE)
  states <- unique(data$state)
  colnames <- c('id', 'inx', 'state', 
                'water_dist', 'evi', "was_burned",
                'pred', 'pred_lo', 'pred_hi')
  
  # set up data frame
  LE=length(evi); LS=length(states); LC=length(colnames);
  pdata <- matrix(0, 
                    nrow=LE*LS,
                    ncol=LC) %>% as.data.frame()
  names(pdata) <- colnames
  
  # fill pdata
  pdata$id <- as.character(id)
  pdata$inx <- 1:nrow(pdata)
  pdata$water_dist <- water_dist
  pdata$evi <- rep(evi, by=LS)
  pdata$state <- rep(states, each=LE)
  
  # predicted fit
  fit <- predict(hsf, pdata, se.fit=TRUE, interval="confidence")
  pdata$pred <- exp(fit$fit)
  pdata$pred_lo <- exp(fit$fit - fit$se.fit)
  pdata$pred_hi <- exp(fit$fit + fit$se.fit)
  
  # return
  return(pdata)
}
pdata <- predictData(id=1, hsf)

# plot to show off the EVI difference between states
ggplot(data=pdata, 
       aes(x=evi, fill=state, color=state) ) + 
  geom_ribbon(mapping=aes(ymin=pred_lo, ymax=pred_hi), 
              linetype='dashed', alpha=0.1) +
  geom_line(mapping=aes(y=pred), linewidth=1)

# plot to show off the responses to different cover types
# ggplot(pdata %>% filter(evi == 0.25), 
#        aes(x=cover_class, fill=state, color=state),
#        position=position_dodge(0.5)) + 
#   geom_segment(aes(y=pred_lo, yend=pred_hi), 
#                lineend="butt", 
#                linewidth=1) +
#   geom_point(aes(y=pred), size=2, color='black')


# plot EVI by cover type
data %>% 
  group_by(cover_class) %>% 
  summarize(evi = mean(evi, na.rm=TRUE), 
            sd  = sd(evi, na.rm=TRUE))
ggplot(data, aes(x=cover_class, y=evi, fill=season)) + 
  geom_violin(linewidth=1, alpha=0.7) + 
  theme(text=element_text(size=16)) + 
  scale_fill_brewer(palette="Dark2", direction=-1)

# ******************************************************************************
#                               SAVE RESULTS
# ******************************************************************************


# I think that hsf3 is the best choice, despite having higher AIC than the 
#   model with a state:cover_class interaction. It's less complicated and 
#   honestly the significance values aren't wowing me.
hsf <- hsf3
save(hsf, file=here(outdir, 'habitat_selection', 'hsf.rdata'))


