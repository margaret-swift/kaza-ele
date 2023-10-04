# hmm.R
# Created 02 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/hmm.R')
source(here('02_scripts','utilities.R'))
pacman::p_load(momentuHMM, rgdal)
setDataPaths('elephant')
load(here(procpath, 'ele.rdata'))

# ******************************************************************************
#                                HMM PREP DATA
# ******************************************************************************

data <- ele.df %>% filter(ID <= 3)

# add step lengths and turning angles again
pdata <- prepData( data,
                   type="UTM",
                   coordNames=c("X", "Y"),
                   covNames=c('SEASON',
                              'SEX',
                              'INX',
                              'DIST'))

hist(pdata$step, breaks=100)

# ******************************************************************************
#                                PARAMETERS
# ******************************************************************************

# label states
nbStates <- 3
stateNames <- c("resting", "foraging", "exploring")

# distributions for observation processes
dist = list(step = "gamma", angle = "vm")

# initial parameters
mu0 <- c(0.1, 0.5, 2)
sd0 <- c(0.1, 0.5, 1)

# plotting hist and params
hist(pdata$step, breaks=100, xlim=c(0, 6000))

# sum(pdata$step == 0, na.rm=TRUE)/nrow(pdata)
zm0 <- c(0.005, 0.001, 0.0005)
angleMean0 <- c(pi,0,0) # initial means (one for each state) 
angleCon0 <- c(0.5, 0.5, 0.85) # initial concentrations (one for each state)

#put it all together
Par0 = list(step = c(mu0, sd0, zm0), 
            angle = c(angleMean0, angleCon0))



# ******************************************************************************
#                             FITTING HMM - SIMPLE
# ******************************************************************************

# formula
form <- formula( ~ TOD + SEASON )

t0 <- proc.time()[['elapsed']]
# fit model
m.roan.simp <- momentuHMM::fitHMM(
  data = pdata.simp %>% filter(SPECIES == "Roan"), 
  nbStates = nbStates, 
  dist = dist, 
  Par0 = Par0,
  estAngleMean = list(angle=TRUE), 
  stateNames = stateNames,
  formula=form)

m.oryx.simp <- fitHMM(
  data = pdata.simp %>% filter(SPECIES == "Oryx"), 
  nbStates = nbStates, 
  dist = dist, 
  Par0 = Par0,
  estAngleMean = list(angle=TRUE), 
  stateNames = stateNames,
  formula=form)
t1 <- proc.time()[['elapsed']]
print(round((t1-t0)/60, 2))

# decoding states and adding to data
roan$STATE <- stateNames[viterbi(m.roan.simp)]
oryx$STATE <- stateNames[viterbi(m.oryx.simp)]
save(m.roan.simp, m.oryx.simp, 
     roan, oryx,
     file=here("03_output", "hmm", "HMMOutputSimple2.RData"))




# ******************************************************************************
#                                PLOTTING TIME!
# ******************************************************************************

# # plot data
plot(roan, compact=T)
plot(oryx, compact=T)
# 
# # show model outputs
# m.roan
# m.oryx
# 
# # plot model outputs
# plot(m.roan.simp, plotCI=TRUE)
plot(m.oryx, plotCI=TRUE)
# 
# plotStates(m, animals="elk-115")
# 
plotStationary(m.roan, plotCI=TRUE)
# plotStationary(m.oryx, plotCI=TRUE)
# 

rbind(roan, oryx) %>%
  group_by(STATE) %>%
  summarize(step=mean(DIST, na.rm=TRUE),
            angle=mean(angle, na.rm=TRUE),
            n=n(),
            p=n/nrow(.))
m %>%
  ggplot() +
  geom_histogram(aes(x=angle)) +
  facet_wrap(~STATE)

