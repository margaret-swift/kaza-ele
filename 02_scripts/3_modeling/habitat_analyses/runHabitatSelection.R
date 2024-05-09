# runHabitatSelection.R
# Created 13 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# amt vignette: 
#   https://cran.r-project.org/web/packages/amt/vignettes/p4_SSF.html
# interpreting SSF results: 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2656.13441


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/habitat_analyses/runHabitatSelection.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(survival)
setDataPaths('geographic')
load(procpath('geographic.rdata'))
load(procpath('stepSelectionParamsHSF.rdata'))
lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))



# ******************************************************************************
#                             RUN CONDITIONAL HSF GLM
# ******************************************************************************

# CONDITIONAL HSF
data <- ssf.df %>% filter(year(t1_) == 2017, month(t1_) %in% 5:9) %>% 
  left_join(lands.meta[,c('category', 'class')], by="category")

# functions
hsfFunc <- function(f, d, type) {
  d$type <- d[,type][[1]]
  clogit(f, 
         data=d, 
         weights=w, 
         method="approximate")
}
hsf1 <- function(d, type) {
  f <- formula(case_ ~ evi + waterdist + type + strata(id))
  hsfFunc(f, d, type)
}
hsf2 <- function(d, type) {
  f <- formula(case_ ~ evi*state + waterdist + type + strata(id))
  hsfFunc(f, d, type)
}
hsf3 <- function(d, type) {
  f <- formula(case_ ~ evi*state + waterdist*state + type + strata(id))
  hsfFunc(f, d, type)
}
hsf4 <- function(d, type) {
  f <- formula(case_ ~ evi*state + waterdist*state + type*state + strata(id))
  hsfFunc(f, d, type)
}

## HSF 1
hsf1.cov <- hsf1(data, 'cover_class')
hsf1.veg <- hsf1(data, 'veg_class')
hsf1.class<-hsf1(data, 'class')

## HSF 2
hsf2.cov <- hsf2(data, 'cover_class')
hsf2.veg <- hsf2(data, 'veg_class')
hsf2.class<-hsf2(data, 'class')

## HSF 3
hsf3.cov <- hsf3(data, 'cover_class')
hsf3.veg <- hsf3(data, 'veg_class')
hsf3.class<-hsf3(data, 'class')

## HSF 4
hsf4.cov <- hsf4(data, 'cover_class')
hsf4.veg <- hsf4(data, 'veg_class')
hsf4.class<-hsf4(data, 'class')


AIC(hsf1.cov, hsf1.veg, hsf1.class,
    hsf2.cov, hsf2.veg, hsf2.class,
    hsf3.cov, hsf3.veg, hsf3.class,
    hsf4.cov, hsf4.veg, hsf4.class
) %>% as.data.frame() %>% arrange(AIC)


## comparing coefficients for different models
getCoefs <- function(hsf) {
  exp(coefficients(hsf)[c('evi', 'waterdist', 'stateforaging', 'stateresting', 'evi:stateresting', 'evi:stateforaging')])
  }
data.frame( cover=getCoefs(hsf2.cov),
            veg = getCoefs(hsf2.veg),
            class=getCoefs(hsf2.class))
summary(hsf2.cov)
summary(hsf2.veg)
summary(hsf2.class)


plotExpCoefs <- function(hsf) {
  cfs <- summary(hsf)$conf.int %>% 
    as.data.frame() %>% 
    rownames_to_column(var='variable')
  names(cfs) <- tolower(gsub(' \\.', '_', 
                             gsub('\\(|\\)', '', names(cfs))))
  
  
  # base state is exploring
  # base landcover is 'bare/water'

  # if you are foraging
  fadj <- cfs$expcoef[cfs$variable=='stateforaging']
  radj <- cfs$expcoef[cfs$variable=='stateresting']
  landcov <- cfs %>% 
    filter(grepl('type', variable)) %>% 
    mutate(variable = gsub('type', '', variable),
           fadjcoef = expcoef * fadj,
           radjcoef = expcoef * radj)
           
  cfs %>% 
    filter(grepl('type', variable)) %>% 
    mutate(variable = gsub('type', '', variable)) %>% 
    ggplot(aes(x=variable, xend=variable)) + 
      geom_hline(yintercept=1, color='darkgray', linetype='dashed') +
      geom_segment(aes(y=lower_95, yend=upper_95)) + 
      geom_point(aes(y=expcoef), size=2) + 
      ylab('difference from bare ground') + xlab('landcover type') +
      theme(text=element_text(size=20))
}


# ******************************************************************************
#                               SAVE RESULTS
# ******************************************************************************


# I think that hsf3 is the best choice, despite having higher AIC than the 
#   model with a state:cover_class interaction. It's less complicated and 
#   honestly the significance values aren't wowing me.
hsf <- hsf3
save(hsf, file=here(outdir, 'habitat_selection', 'hsf.rdata'))


