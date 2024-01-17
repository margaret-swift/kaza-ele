# determineActivityStates.R
# Created 02 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# momentuhmm vignette: https://cran.r-project.org/web/packages/momentuHMM/vignettes/momentuHMM.pdf

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
pacman::p_load(momentuHMM, rgdal, here, ConnMatTools)
i_am('02_scripts/3_modeling/determineActivityStates.R')
source(here('02_scripts', 'utilities.R'))
setDataPaths('elephant')
load(procpath('ele.rdata'))
load(here(outdir, 'hmmSeasonF.rdata'))
setDataPaths('precipitation')
load(procpath('precipitation.rdata'))

# ******************************************************************************
#                                HMM PREP DATA
# ******************************************************************************

# choose only certain IDs for a smaller run
ids <- unique(ele.df$ANIMAL_ID)
data <- ele.df %>% nog() %>% 
  dplyr::filter(SEX == "F", 
                ANIMAL_ID %in% ids[1:12],
                FIXRATE == "1 hour") %>%
  mutate(SEASON = factor(SEASON),
         ID=BURST) %>% 
  dplyr::select(INX, ID, X, Y, 
                SEASON, SEX, DIST, DATE.TIME)

# Crawl Wrap: Fill in missing data steps
crwOut <-crawlWrap(obsData = data,
                   timeStep = "hour",
                   Time.name = "DATE.TIME",
                   coord=c('X', 'Y'),
                   ncore=3
                   )
crwData <- crwOut$crwPredict 
crwData$date <- as.POSIXct(crwData$DATE.TIME)

# set SEASON based on first rainfall from precipitation.rdata
crwData$ISWET <- TRUE
crwData$YEAR <- year(crwData$date)
dt <- 121
inx <- which( yday(crwData$date) > dt )
inx.rain <- match(crwData$YEAR+1, wet.start$RAINYEAR)
iswet <- crwData$date > wet.start$DATE[inx.rain]
crwData$ISWET[inx] <- iswet[inx]
crwData$SEASON = ifelse(crwData$ISWET, "WET", "DRY")
crwData <- crwData[,c('ID', 'X', 'Y', 'SEASON', 'date')]

# Prep Data: add step lengths and turning angles again
tic()
pdata <- prepData( crwOut,
                   type="UTM",
                   coordNames=c("X", "Y"),
                   covNames=c('SEASON')
                   )
toc()

# ******************************************************************************
#                                PARAMETERS
# ******************************************************************************

# label states
nbStates <- 3
stateNames <- c("resting", "foraging", "correlated walk")

# distributions for observation processes
dist = list(step = "gamma", angle = "wrpcauchy")

# STEP LENGTH
mu0.s <- c(10, 400, 2000)
sd0.s <- c(20, 300, 1000)
zm0.s <- c(0.05, 0.001, 0.0005)
plotParams(mu0.s, sd0.s, "gamma")

# TURNING ANGLE
mu0.a <- c(pi, pi, 0) # initial means (one for each state)
kappa0.a <- c(0.1, 0.5, 2) # initial concentrations (one for each state)
plotParams(mu0.a, kappa0.a, "wrpcauchy")

#put it all together
Par0 = list(step = c(mu0.s, sd0.s, zm0.s), 
            angle = c(mu0.a, kappa0.a))

# formula
form <- formula( ~ SEASON )


# ******************************************************************************
#                             FITTING HMM - SIMPLE
# ******************************************************************************

tic()
# fit model
hmm <- momentuHMM::fitHMM(
  data = pdata,
  nbStates = nbStates, 
  dist = dist, 
  Par0 = Par0,
  estAngleMean = list(angle=TRUE), 
  # circularAngleMean=list(angle=TRUE),
  stationary = FALSE,
  stateNames = stateNames,
  modelName = "female_with_seasons",
  formula=form)
toc()


checkPar0(pdata, nbStates, dist, Par0)

# # save MALE data
# pdata$STATE <- data$STATE <- stateNames[viterbi(hmm)]
# pdata.m <- pdata
# data.m <- data
# hmm.m <- hmm
# save(hmm.m, pdata.m, data.m, file=here(outdir, 'hmmSeasonM.rdata'))

# save FEMALE data
pdata$STATE <- data$STATE <- stateNames[viterbi(hmm)]
pdata.f <- pdata
data.f <- data
hmm.f <- hmm
save(hmm.f, pdata.f, data.f, file=here(outdir, 'hmmSeasonF.rdata'))


# ******************************************************************************
#                                PLOTTING TIME!
# ******************************************************************************

# # plot data
plot(data, compact=T)

# show model outputs
hmm <- hmm.f

# plot model outputs
plot(hmm, plotCI=TRUE)

# plot stationary data
plotStationary(hmm, plotCI=TRUE)

pdata %>%
  ggplot() +
  geom_histogram(aes(x=angle)) +
  facet_wrap(~STATE)


