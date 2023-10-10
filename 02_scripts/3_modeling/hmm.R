# hmm.R
# Created 02 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
pacman::p_load(momentuHMM, rgdal, here)
i_am('02_scripts/3_modeling/hmm.R')
source(here('02_scripts','utilities.R'))
setDataPaths('elephant')
load(here(procpath, 'ele.rdata'))
load(here(outdir, 'hmm.rdata'))

# ******************************************************************************
#                                HMM PREP DATA
# ******************************************************************************

data <- ele %>% nog() %>% 
  filter(ID == 1,
         FIXRATE == "1 hour", 
         YEAR == 2016) %>%
  dplyr::select(INX, SEASON, ANIMAL_ID, SEX, DIST,
                DATE.TIME, BURST, X, Y) %>% 
  rename(ID=BURST)

data$DIFF <- as.n(abs(difftime(lead(data$DATE.TIME), data$DATE.TIME)))
data <- data %>% 
  mutate(X = ifelse(DIFF>100, NA, X),
         Y = ifelse(DIFF>100, NA, Y))

# add step lengths and turning angles again
pdata <- prepData( data,
                   type="UTM",
                   coordNames=c("X", "Y"),
                   covNames=c('SEASON',
                              'ANIMAL_ID',
                              'BURST',
                              'SEX',
                              'INX',
                              'DIST'))
which(pdata$step > max(data$DIST, na.rm=TRUE))
pdata[3155:3157,]

# ******************************************************************************
#                                PARAMETERS
# ******************************************************************************

# label states
nbStates <- 3
stateNames <- c("resting", "foraging", "correlated walk")

# distributions for observation processes
dist = list(step = "gamma", angle = "vm")

# initial parameters
mu0 <- c(10, 200, 1500)
sd0 <- c(5, 300, 1000)
# mu0 <- c(200, 1500)
# sd0 <- c(200, 1200)

# plotting hist and params
hist(pdata$step, breaks=200)
abline(v=mu0, col='red')
abline(v=mu0+sd0, col='pink')
abline(v=mu0-sd0, col='pink')

zm0 <- c(0.005, 0.001, 0.0005)
angleMean0 <- c(pi,0,0) # initial means (one for each state)
angleCon0 <- c(0.5, 0.5, 0.85) # initial concentrations (one for each state)
# zm0 <- c(0.001, 0.0005)
# angleMean0 <- c(pi,0) # initial means (one for each state) 
# angleCon0 <- c(0.5, 0.85) # initial concentrations (one for each state)

#put it all together
Par0 = list(step = c(mu0, sd0),# zm0), 
            angle = c(angleMean0, angleCon0))



# ******************************************************************************
#                             FITTING HMM - SIMPLE
# ******************************************************************************

# formula
form <- formula( ~ 1 )

tic()
# fit model
hmm <- momentuHMM::fitHMM(
  data = pdata, # pdata %>% filter(SEX=="F")
  nbStates = nbStates, 
  dist = dist, 
  Par0 = Par0,
  estAngleMean = list(angle=TRUE), 
  stateNames = stateNames,
  formula=form)
toc()

# decoding states and adding to data
pdata$STATE <- stateNames[viterbi(hmm)]
data$STATE <- stateNames[viterbi(hmm)]

save(hmm, file=here(outdir, 'hmmShort.rdata'))


# ******************************************************************************
#                                PLOTTING TIME!
# ******************************************************************************

# # plot data
plot(data, compact=T)

# show model outputs
hmm

# plot model outputs
plot(hmm, plotCI=TRUE)
points(data$X[data$STATE=="resting"],
     data$Y[data$STATE=="resting"],
     col='orange', cex=0.5, pch=19)

plotStationary(m.roan, plotCI=TRUE)

pdata %>%
  ggplot() +
  geom_histogram(aes(x=angle)) +
  facet_wrap(~STATE)

