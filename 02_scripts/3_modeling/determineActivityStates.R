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

# choose only females and 1-hour fix rates
data <- ele.df %>% nog() %>% 
  dplyr::filter(SEX == "F", 
                FIXRATE == "1 hour") %>%
  mutate(SEASON=factor(SZN_4),
         ID=BURST) %>% 
  dplyr::select(INX, ID, X, Y, 
                SEASON, SEX, DIST, DATE.TIME)

# Crawl Wrap: Fill in missing data steps
tic()
crwOut <-crawlWrap(obsData = data,
                   timeStep = "hour",
                   Time.name = "DATE.TIME",
                   coord=c('X', 'Y'),
                   ncore=3
                   )
crwOut$crwPredict$date <- as.POSIXct(crwOut$crwPredict $DATE.TIME)
toc()

# Prep Data: add step lengths and turning angles again
tic()
# use top ten eles with most data
crw <- crwOut
idlist <- names(sort(table(crw$crwPredict$ID))[23:33]) %>% as.numeric()
crw$crwPredict <- crw$crwPredict[crw$crwPredict$ID %in% idlist,]
crw$crwFits <- crw$crwFits[names(crw$crwFits) %in% idlist]

# run pdata
pdata <- prepData( crw,
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
sd0.s <- c(10, 300, 1500)
zm0.s <- c(0.5, 0.01, 0.005)
plotParams(mu0.s, sd0.s, "gamma")

# TURNING ANGLE
mu0.a <- c(0, pi, 0) # initial means (one for each state)
kappa0.a <- c(0.001, 1, 1) # initial concentrations (one for each state)
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

colors <- c('#408F67', '#61AAE1', '#D99B35')
pdata %>%
  ggplot() +
  geom_histogram(aes(x=angle, fill=STATE), bins=30) +
  facet_wrap(~STATE) + 
  scale_fill_manual(values=colors) + plot.theme + 
  xlab("angle (radians)") + ylab('number of points')

pdata %>%
  ggplot() +
  geom_histogram(aes(x=step, fill=STATE), bins=100) +
  facet_wrap(~STATE) + 
  scale_fill_manual(values=colors) +
  plot.theme + 
  xlab("step size (m)") + ylab('number of points')

pdata %>%

pangle / pstep + plot_layout(guides="collect")


