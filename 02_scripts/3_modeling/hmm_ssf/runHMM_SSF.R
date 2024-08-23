# runHMM_SSF.R
# Created 02 June 2024
# Margaret Swift <margaret.swift@cornell.edu>

# HMM SSF package info
# devtools::install_github("NJKlappstein/hmmSSF")
# hmmSSF vignette: 
# https://github.com/NJKlappstein/hmmSSF/blob/main/vignettes/hmmSSF_introduction.pdf

#                     === Word of caution from the authors ===
# This model is fairly complex compared to standard step selection functions, and 
#   requires more care to avoid numerical problems. Indeed, the number of parameters 
#   is higher, and the problem of inference requires estimating the dynamics of an 
#   unobserved state process (used as proxy for the behavioural state of the animal). 
# 
# In our experience, this model complexity often leads to numerical instability, 
#   i.e., failure of the fitting function to converge to the best parameter estimates. 
#   This is a difficult general problem, and here are some possible directions to 
#   explore to help with it:
#   
#   - increase the number of control locations
#   - try different sets of starting values for the parameters
#   - use a simpler model formulation. The more complex the model, the more 
#       difficult it is to fit. It is usually good to start from the simplest 
#       model, i.e., a 2-state model with no covariates on the transition 
#       probabilities, and build up complexity if it seems possible.
#   
# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

# packages and location
pacman::p_load(here, hmmSSF)
i_am('02_scripts/3_modeling/hmm_ssf/runHMM_SSF.R')

# utilities
source(here('02_scripts', 'utilities.R'))

# loading data
setOutPath('hmm')
load(outpath('hmm_ssf.rdata'))
data <- hmm.ssf %>% 
  mutate(ID=id, #wtf i can't win
         step = step / 1000) #tutorial is in km not m
rm(hmm.ssf)
gc()

# remove EVI NA's and replace with the mean
emean <- mean(data$evi)
data1 <- data %>% mutate(evi = ifelse(is.na(evi), emean, evi))
rm(data)
gc()

# ******************************************************************************
#                       ==== STARTING PARAMETER VALUES ====
# ******************************************************************************
#                     

# In hmmSSF, the model is fitted through numerical likelihood optimization (nlm()), 
#   and starting values need to be chosen for the model parameters. Poorly chosen 
#   values can lead to numerical problems (e.g., failure to converge). 
# If starting values are close to the "true" parameters,it will be easier for the 
#   optimiser to converge. We choose plausible values given the data. 
# Here, we use the fact that the selection parameters for step and log(step) are 
#   linked to the mean and sd of a gamma distribution. We can choose some 
#   hypothesised mean step length for each state by looking at the empirical 
#   distribution, and transform them using these formulas:

betaL = function(muL, sigL) { - ( muL / (sigL^2) ) }
betaLogL = function(muL, sigL) { ( (muL^2) / (sigL^2) ) - 2 }

# We want state 1 to capture slow movement and state 2 to capture fast movement:
s1mu <- s1sd <- 0.1
s2mu <- s2sd <- 0.7
s3mu <- s3sd <- 2
betaLParams_3    <- c(betaL(s1mu, s1sd), 
                    betaL(s2mu, s2sd),
                    betaL(s3mu, s3sd))
betaLogLParams_3 <- c(betaLogL(s1mu, s1sd), 
                    betaLogL(s2mu, s2sd),
                    betaLogL(s3mu, s3sd))
betaLParams_2 <- betaLParams_3[c(1,3)]
betaLogLParams_2 <- betaLogLParams_3[c(1,3)]

# Similarly, the selection parameter for cos(angle) is the concentration of a 
#   von Mises distribution. The larger the concentration, the more directed the 
#   movement. We choose a starting value of 0.2 in state 1 and 5 in state 2.
kappa_3 <- c(0.2, 1, 5)
kappa_2 <- c(0.2, 5)

# The parameters should be passed as a matrix, with one row for each covariate 
#   in the formula (here, 4), and one column for each state (here, 2).
evi.params_3 <- c(2, 4, -1)
evi.params_2 <- c(3, -1)

# put em all together
dat_3 <- c(betaLParams_3,
           betaLogLParams_3,
           kappa_3, # concentration
           evi.params_3) # evi... what about factors (veg class?)
dat_2 <- c(betaLParams_2,
           betaLogLParams_2,
           kappa_2,
           evi.params_2)

# ******************************************************************************
#                ==== SET NUMBER OF STATES & SAVE VALUES ====
# ******************************************************************************

# There is no general method to select the best number of states, and we recommend 
#   leaning towards 2 or 3 states in most applications to avoid numerical issues.
n_states <- 2
if (n_states == 3) {dat=dat_3
} else {dat=dat_2}

# initial values for SSF parameters
ssf_par0 <- matrix(dat, ncol=n_states, byrow=TRUE)

# The formula should include movement variables as well as covariates (@avgar2016). 
ssf_formula <- ~ step + log_sl + cos_ta + evi

# Getting an ERROR: function cannot be evaluated at initial parameters
#  removed EVI for the moment but this might help? 
#  https://stackoverflow.com/questions/69850540/error-in-optim-r-cannot-be-evaluated-at-initial-parameters

data1$animal_id <- gsub('\\d{3}$', '', data1$id)
data2 <- data1 %>% filter(animal_id %in% 1:10)

# Finally, we can pass the above to fit the model.
mod <- hmmSSF(ssf_formula = ssf_formula,
              n_states = n_states,
              data = data2,
              ssf_par0 = ssf_par0)

# STS
save(mod, file=outpath('hmm_ssf_mod.rdata'))













# # To visualise step lengths, we subset the data to observed steps (i.e., removing 
# #   control steps).
# ggplot(subset(data, obs == 1), aes(x = step)) +
#   geom_histogram(col = "white", fill = "grey", bins = 30) +
#   # geom_line(data=s1dist, mapping=aes(x=x, y=y)) +
#   # geom_line(data=s2dist, mapping=aes(x=x, y=y)) +
#   geom_vline(aes(xintercept = s1mu), linewidth=1, linetype='dashed', color='purple') +
#   geom_vline(aes(xintercept = s2mu), linewidth=1, linetype='dashed', color='blue') +
#   geom_vline(aes(xintercept = s3mu), linewidth=1, linetype='dashed', color='green') +
#   big.theme
#   
#   
#   
#   # createGamma <- function(mu, sd) {
#   x <- seq(0, 6000, length.out=100)
#   shape = mu^2 / sd^2
#   scale = sd^2 / mu
#   y <- dgamma(x, shape=shape, scale=scale)
#   return(data.frame(x,y))
# }
# s1dist <- createGamma(s1mu, s1sd)
# s2dist <- createGamma(s2mu, s2sd)