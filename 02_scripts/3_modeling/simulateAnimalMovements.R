# simulateAnimalMovements.R
# Created 16 Oct 2023
# Margaret Swift <mes114@duke.edu>

# Marshall BM and Duthie AB. abmAnimalMovement: An R package for simulating 
#   animal movement using an agent-based model [version 1; peer review: 2 
#   approved with reservations]. F1000Research 2022, 11:1182 
#   (https://doi.org/10.12688/f1000research.124810.1)

# devtools::install_github("BenMMarshall/abmAnimalMovement")
# GitHub page: https://github.com/BenMMarshall/abmAnimalMovement
# Vignettes: https://joshcullen.github.io/bayesmove/articles/index.html

# "The model simulates animal locations over a given period of time at discrete 
#  time steps. At each time step, the agent (i.e., simulated animal) is presented
#  with a range of movement options in the form of new locations [Figure 2], 
#  and will choose from amongst these (i.e., sum-based model as opposed to 
#  facing a series of sequential binary decisions, see Ref. 2 for an example 
#  of the latter). The possible movement options, and how the new location 
#  is selected, are influenced by several factors: behavioural state of the 
#  animal, environmental quality, and proximity to points of attraction/avoidance. 
#  
#  By simulating these drivers of animal movement, we hope to capture aspects 
#  of the:
#  - internal state       (i.e., motivation to move via behaviour); 
#  - motion capacity      (i.e., the individual’s varying ability to move); 
#  - navigation capacity  (i.e., ability to plan ahead beyond the immediate 
#                                movement distance); and
#  - external factors     (i.e., the landscape that steers and limits movement), 
#  all of which are defined as key components of animal movement."



# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

pacman::p_load(here,
       furrr,
       dplyr,
       nngeo,
       reshape2,
       ggplot2,
       ggforce,
       ggtext,
       ggridges,
       patchwork,
       raster)
       # NLMR)
source(here('02_scripts','utilities.R'))
setDataPaths('elephant')
load(procpath('ele.rdata'))

setDataPaths('geographic')
load(procpath('geographic.rdata'))
load(procpath('landsmats.rdata'))

# install abmFences
roxygen2::roxygenize("~/R_Packages/abmFences")
# devtools::install_github("margaret-swift/abmFences")
# library(abmFences)

# set flag for which sex to simulate
sex <- "F" #"M"

# sex-specific parameters
# Real fence encounter rate data from Naidoo et al 2022, at 1km threshold
# for fences, rivers, and roads
if (sex == "F") {
  load(here(outdir, 'hmmLongF.rdata'))
  mle <- hmm.f$mle
  perm <- c(0, 0.101, 0.153)
} else { 
  load(here(outdir, 'hmmLongM.rdata')) 
  mle <- hmm.m$mle
  perm <- c(0.035, 0.145, 0.258)
}

# ******************************************************************************
#                     Get transition matrix and params from HMM
# ******************************************************************************

ELE_behaveMatrix <- matrix(nrow=3, ncol=3)
for (i in 1:3) {
  # get transition matrix data from HMM MLE
  d2 <- exp(mle$beta[,(i*2)-1])
  d3 <- exp(mle$beta[,(i*2)])
  #calculate probs
  p2 <- d2 / ( 1 + d2 + d3 )
  p3 <- d3 / ( 1 + d2 + d3 )
  p1 <- 1 - ( p2 + p3 )
  #save
  ELE_behaveMatrix[i,] <- c( p1, p2, p3 )
}

# step length (CPP CODE ASSUMES MINUTE SCALE)
ELE_k_step <- mle$step['mean',] / 60
ELE_s_step <- mle$step['sd',] / 60

# turning angle
ELE_mu_angle <- mle$angle['mean',]
ELE_k_angle  <- mle$angle['concentration',]


# ******************************************************************************
#                         Set up raster and barrier layers
# ******************************************************************************

ELE_shelter<- shelter.mat
ELE_forage <- forage.mat
ELE_move   <- move.mat

# barriers -- skip roads for now
barriers <- list(fences, rivers)#, roads)
barriers <- sapply(barriers, function(e) st_transform(e, crs=32735))
barriers_data <- generateBarriers(barriers, perm[1:2])


# ******************************************************************************
#                         Remaining simulation parameters
# ******************************************************************************

# shelters
ELE_shelterLocs <- data.frame(
  "x" = c(90, 95, 110),
  "y" = c(90, 75, 90))
ELE_shelterSize <- 5

# destination params
ELE_destinationRange <- c(3, 120)
ELE_destinationDirection <- c(0, 0.01)
ELE_destinationTransformation <- 2
ELE_destinationModifier <- 2

# scale of raster resolution
ELE_rescale <- mat.res[1] #should be around 10

# avoiding
# avoidTrans 0 - no transformation applied to the distance to avoidance 
#                 points weighting, 
#            1 - distance to avoidance points weighing is square-rooted, 
#            2 - distance to avoidance points weighting is squared
# avoidMod  A coefficient to be applied to the avoidance points weighting.
ELE_avoidLocs <- data.frame( "x" = 5, "y" = 5 )
ELE_avoidTransformation <- 2
ELE_avoidModifier <- 4

# additional cycles
ELE_rest_Cycle <- c(0.12, 0, 24, 24)
c0 <- c(0.075, 0, 24* (365/2), 24* 365) # seasonal
ELE_additional_Cycles <- rbind(c0)

# choose start location and options
start <- c(26, -20)
timesteps <- 5000 #24*60*31
des_options=10; options=12;

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                         SIMULATE ELEPHANT MOVEMENTS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

set.seed(1001)

## Run this as a chunk to generate random maps of elephant movements
seed <- floor(runif(1, 0, 1) * 1e5)
runSim(seed, barriers, barriers_data, perm, colorby='inx')


















# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






# row <- 200; col <- 200;
# landscapeLayersList <- list(
#   "shelter" = matrix(runif(row*col, 0, 1), nrow = row, ncol = col),
#   "forage"  = matrix(runif(row*col, 0, 1), nrow = row, ncol = col),
#   "movement"= matrix(runif(row*col, 0, 1), nrow = row, ncol = col))
# 
# ELE_shelter <- landscapeLayersList$shelter
# ELE_forage <- landscapeLayersList$forage
# ELE_move <- landscapeLayersList$movement