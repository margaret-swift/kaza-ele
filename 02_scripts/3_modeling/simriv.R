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
i_am('02_scripts/3_modeling/simriv.R')
source(here('02_scripts','utilities.R'))
p_load("SiMRiv", #v1.0.6
       "raster") 
load(here(outdir, 'hmmShort.rdata'))

# ******************************************************************************
#                           MESSING AROUND WITH SIMRIV
# ******************************************************************************

# define a species with a single-state movement type
# characterized by a random walk
SL = 15
# rand.walker <- species(state.CRW(0.98) + SL)
# sim.rw <- simulate(rand.walker, 10000)
# plot(sim.rw, type = "l", asp = 1, main = "Random walk")

# Levy walks
# get estimates for all transitions
vals <- 1:3
probs.df <- data.frame(matrix(nrow=1, ncol=9))
names(probs.df) <- paste(rep(vals, each=3), vals, sep=" -> ")
mat <- hmm$mle$beta
for (i in vals) {
  types <- paste(i, vals[vals != i], sep=" -> ")
  
  d2 <- exp(mat[,types[1]])
  d3 <- exp(mat[,types[2]])
  sum <- d2 + d3
  
  p2 <- d2 / ( 1 + sum )
  p3 <- d3 / ( 1 + sum )
  p1 <- 1 - ( p2 + p3 )
  
  probs <- data.frame( p1, p2, p3 )
  names(probs) <- c(paste(i, i, sep=" -> "), types)
  probs.df[,names(probs)] <- probs
}
tm <- matrix(probs.df, nrow=3)
tm <- transitionMatrix(0.005, 0.01, #state 1>2 and 1>3
                       0.05,  0.02, #state 2>1 and 2>3
                       0.05,  0.005)#state 3>1 and 3>2
levy.walker <- species(
  state.RW() + state.CRW(0.98) + state.Resting(), 
  trans = tm) + SL
sim.lw <- simulate(levy.walker, 1e5)
plot(sim.lw, type='l', asp=1, main='levy walker')
