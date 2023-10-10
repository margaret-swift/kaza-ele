# simriv.R
# Created 04 Oct 2023
# Margaret Swift <mes114@duke.edu>

# Quaglietta, L. & Porto, M. SiMRiv: an R package for mechanistic simulation 
#     of individual, spatially-explicit multistate movements in rivers, 
#     heterogeneous and homogeneous spaces incorporating landscape bias. 
#     Mov Ecol 7, 11 (2019).
# GitHub page: https://github.com/miguel-porto/SiMRiv
# Vignette: https://cran.r-project.org/web/packages/SiMRiv/vignettes/SiMRiv.pdf

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

library(here)
<<<<<<< HEAD
i_am('02_scripts/3_modeling/simriv.R')
source(here('02_scripts','utilities.R'))
p_load("SiMRiv", #v1.0.6
       "raster") 
load(here(outdir, 'hmmShort.rdata'))
=======
i_am(here('02_scripts', 'modeling', 'simriv.R'))
source(here('02_scripts','utilities.R'))
p_load("SiMRiv", #v1.0.6
       "raster") 
>>>>>>> 153847358d2badc07d7c2857d9e2e749a2ce53ca

# ******************************************************************************
#                           MESSING AROUND WITH SIMRIV
# ******************************************************************************

# define a species with a single-state movement type
# characterized by a random walk
SL = 15
<<<<<<< HEAD
# rand.walker <- species(state.CRW(0.98) + SL)
# sim.rw <- simulate(rand.walker, 10000)
# plot(sim.rw, type = "l", asp = 1, main = "Random walk")


# Levy walks
tm <- transitionMatrix( as.n( hmm$mle$beta ) )
=======
rand.walker <- species(state.CRW(0.98) + SL)
sim.rw <- simulate(rand.walker, 10000)
plot(sim.rw, type = "l", asp = 1, main = "Random walk")


# Levy walks
>>>>>>> 153847358d2badc07d7c2857d9e2e749a2ce53ca
tm <- transitionMatrix(0.005, 0.01, #state 1>2 and 1>3
                       0.05,  0.02, #state 2>1 and 2>3
                       0.05,  0.005)#state 3>1 and 3>2
levy.walker <- species(
  state.RW() + state.CRW(0.98) + state.Resting(), 
  trans = tm) + SL
sim.lw <- simulate(levy.walker, 1e5)
plot(sim.lw, type='l', asp=1, main='levy walker')
