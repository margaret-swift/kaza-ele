# runHabitatSelection.R
# Created 13 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# amt vignette: 
#   https://cran.r-project.org/web/packages/amt/vignettes/p4_SSF.html
# interpreting SSF results: 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2656.13441




## NOTES
## -- State foraging, exploring, resting values = baseline use vs. available, 
##    I think this just has to do with the prevalance of those states in the 
##    dataset. 

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/habitat_analyses/runHabitatSelection.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(survival)
setDataPaths('geographic')
load(procpath('geographic.rdata'))
load(procpath('stepSelectionParamsHSF.rdata'))

# join landscape metadata class names to ssf data frame
lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv')) %>% 
  mutate(category = as.factor(category))
lands.meta$newcat <- c('bare', 'bare', 'built',
                       'closed', 'closed', 'closed', 'closed',
                       'cropland', 'open bush', 'open herb', 
                       'open wetland', 'open wetland', 
                       'open bush', 
                       'sparse bush', 'sparse bush', 'sparse bush', 
                       'water', 'water'
                       )

data <- ssf %>% 
  left_join(lands.meta[,c('category', 'class', 'newcat')], by="category")
  # mutate(id = as.factor(id)) 
  



# ******************************************************************************
#                             RUN CONDITIONAL HSF GLM
# ******************************************************************************

# We are using the clogit() function from the "survival" package, where case 1
# is the step actually taken and case 0 is steps not taken. 

# CONDITIONAL HSF
hsfFunc <- function(f, data, type) {
  data$type <- data[,type][[1]]
  clogit(f, 
         data=data, 
         weights=w, 
         method="approximate")
}
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
# MODEL SELECTION
modelSelection <- function(d, type) {
  f0 <- formula(case_ ~ evi + water_dist + was_burned + strata(inx))
  f1 <- formula(case_ ~ evi*state + water_dist*season + strata(inx))
  f2 <- formula(case_ ~ evi*state + was_burned + water_dist+ strata(inx))
  f3 <- formula(case_ ~ evi*state + water_dist + was_burned + type + strata(inx))
  hsf_ls <- list( hsfFunc(f0, d, type),
                  hsfFunc(f1, d, type),
                  hsfFunc(f2, d, type),
                  hsfFunc(f3, d, type))
  return(hsf_ls)
}

## HSF selection
hsf.cov <- modelSelection(data, 'newcat')
# hsf.veg <- modelSelection(data, 'veg_class')
hsf.class<-modelSelection(data, 'class')

# GET AIC AND ARRANGE NICELY
compareModels <- function(mlist, name) {
  calls <- sapply(mlist, function(e) as.character(e$formula[3]))
  df <- data.frame(model=paste(name, calls, sep=" : "),
                   AIC=sapply(mlist, AIC), 
                   LRT=sapply(mlist, function(e) summary(e)$logtest['test']),
                   LRTp=sapply(mlist, function(e) summary(e)$logtest['pvalue']),
                   rsq= sapply(mlist, function(e) summary(e)$rsq[[1]]))
  return(df)
}
# LISTING AIC
aic.df <- rbind(compareModels(hsf.cov, "cover_class"), 
                compareModels(hsf.class, "class")) %>% 
  arrange(AIC)
aic.ref <- aic.df$AIC[1]
aic.df$dAIC <- round(aic.df$AIC - aic.ref)
aic.df

## comparing coefficients for different models
for (i in 1:length(hsf.cov)) {
  spacer<-"\n============\n"
  message(spacer, i, spacer)
  x1 <- hsf.cov[[i]]
  c1 <- getCoefs(x1)
  daic1 <- round(AIC(x1)-aic.ref)
  x2 <- hsf.class[[i]]
  c2 <- getCoefs(x2)
  daic2 <- round(AIC(x2)-aic.ref)
  if (nrow(c1) != nrow(c2)) {
    message('cover type - dAIC = ', daic1)
    print(c1)
    message('\nfull landcover class - dAIC = ', daic2)
  } else { message('dAIC = ', daic1)}
  print(c2)
}




# ******************************************************************************
#                     RUN MODELS SEPARATELY FOR EACH STATE
# ******************************************************************************

# answers for the error "loglik converged before variable N"
# https://stat.ethz.ch/pipermail/r-help/2008-September/174201.html

hsfSznState <- function(f, d, st, type='newcat') {
  szns <- c('dry', 'wet')
  mylist <- list('dry'=NULL, 'wet'=NULL)
  d_st <- d[d$state == st,]
  
  for (szn in szns) {
    # run dry data
    d_sz <- d_st[d_st$season == szn,]
    NR <- nrow(d_sz)
    run <- paste(szn, st)
    mod <- NULL
    message('modeling ', run)
    message("  # rows = ", NR)
    if (NR > 10000) {
      mod =  hsfFunc(f, data=d_sz, type=type)
      mylist[[szn]] <- mod
    } else { message('WARNING: not enough data to run ', run) }
  }
  return(mylist)
}
collateMe <- function(mod, st, szn) {
  ncol <- 3
  if (!is.null(mod)) {
    types <- paste0('newcat', unique(data$newcat))
    df <- getCoefs(mod)
    inx.add <- which(!(types %in% rownames(df)))
    df[types[inx.add],] <- rep(NA, ncol)
    df <- df[order(rownames(df)),]
    df <- rownames_to_column(df)
  } else { df = data.frame(matrix(NA, nrow=10, ncol=ncol))}
  names(df) <- c('covar', 'value', 'se', 'signf')
  df <- cbind(state=st, season=szn, df)
  df
}

# set formula and list to hold results
f <- formula(case_ ~ evi + water_dist + newcat + strata(inx))
ids <- 1:4 #unique(data$id)
NL <- length(ids)
summ.mega.ls <- vector(mode='list', length=NL)
names(summ.mega.ls) <- ids

# run model for all together
d <- data[id %in% ids,]
hsf.r <- hsfSznState(f, data, st='resting')
hsf.f <- hsfSznState(f, data, st='foraging')
hsf.e <- hsfSznState(f, data, st='exploring')


summ.all <- rbind(collateMe(hsf.r$dry, 'rest', 'dry'), 
                 collateMe(hsf.f$dry, 'forage', 'dry'), 
                 collateMe(hsf.e$dry, 'explore', 'dry'), 
                 collateMe(hsf.r$wet, 'rest', 'wet'), 
                 collateMe(hsf.f$wet, 'forage', 'wet'), 
                 collateMe(hsf.e$wet, 'explore', 'wet'))

# loop over ids and run model for each state and season.
for (i in 1:NL) {
  id = ids[i];
  d_i <- data[data$id == id,]
  message('RUNNING MODELS FOR ID: ', id)
  hsf.r <- hsfSznState(f, d_i, st='resting')
  hsf.f <- hsfSznState(f, d_i, st='foraging')
  hsf.e <- hsfSznState(f, d_i, st='exploring')
  
  summ.ls <- rbind(collateMe(hsf.r$dry, 'rest', 'dry'), 
                  collateMe(hsf.f$dry, 'forage', 'dry'), 
                  collateMe(hsf.e$dry, 'explore', 'dry'), 
                  collateMe(hsf.r$wet, 'rest', 'wet'), 
                  collateMe(hsf.f$wet, 'forage', 'wet'), 
                  collateMe(hsf.e$wet, 'explore', 'wet'))
  summ.mega.ls[[i]] <- summ.ls
}

# plot cover usage
plotCover <- function(dat) {
  dat$min <- dat$value - dat$se
  dat$max <- dat$value + dat$se
  dat <- dat[grepl('newcat', dat$covar),] %>% 
    mutate(covar = gsub('newcat', '', covar)) %>% 
    filter(max < 5)
  ggplot(dat, aes(color=season, group=season, x=covar)) + 
    geom_linerange(aes(ymin=min, ymax=max),
                   linewidth=2, alpha=0.4,
                   position = position_dodge2(width = 0.5)) +
    geom_point(aes(y=value), size=4, 
               position = position_dodge2(width = 0.5))+
    geom_hline(yintercept=1.0, linetype='dashed', 
               color='darkgray') + 
    facet_wrap(~state, scales='free_y') + 
    scale_color_brewer(palette="Dark2", direction=-1)
}

# for total model
plotCover(summ.all)


i=1
id <- ids[i]
s_i <- summ.mega.ls[[id]]
plotCover(s_i)




plotEstimates <- function(hsf, name) {
  
}

# ******************************************************************************
#                               ANALYZE RESULTS
# ******************************************************************************

# we are going with state x evi and water distance
f <- formula(case_ ~ evi*state + waterdist + strata(id))
hsf <-  clogit(f, 
               data=data, 
               weights=w, 
               method="approximate")

hsf <- hsf.cov[[1]]
summary(hsf)

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
..glm.ratio(hsf)


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


