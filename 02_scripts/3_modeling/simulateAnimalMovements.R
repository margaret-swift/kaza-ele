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
       dplyr,
       ggplot2,
       patchwork)
source(here('02_scripts','utilities.R'))
setDataPaths('elephant', verbose=FALSE)
load(procpath('ele.rdata'))

## Load and or generate raster and barrier layers
setDataPaths('geographic', verbose=FALSE)
load(procpath('geographic.rdata'))
shelt.rast <- terra::rast(procpath('shelterRaster.tif'))
move.rast <- terra::rast(procpath('movementRaster.tif'))
forage.rast <- terra::rast(procpath('forageRaster.tif'))
barriers <- list(fences, rivers)
rasters <- list(shelt.rast, move.rast, forage.rast)
# ELE_landsdata <- generateLandscapeData(rasters, barriers, perms=c(0, 0.101))
# save(ELE_landsdata, file=procpath('simulationInputData.rdata'))
load(procpath('simulationInputData.rdata'))

# build abmAME
# roxygen2::roxygenize("~/R_Packages/abmAME")
# devtools::install_github("margaret-swift/abmAME")
library(abmAME)

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
#                         Remaining simulation parameters
# ******************************************************************************

# shelters
ELE_shelterLocs <- data.frame(
  "x" = c(-155000, -160000, -125000),
  "y" = c(7950000, 7945000, 7945000))
ELE_shelterSize <- 50

# destination params
ELE_destinationRange <- c(3, 120)
ELE_destinationDirection <- c(0, 0.01)
ELE_destinationTransformation <- 2
ELE_destinationModifier <- 2

# scale of raster resolution
ELE_rescale <- ELE_landsdata$landscape_meta[5] #should be around 10

# avoiding
# avoidTrans 0 - no transformation applied to the distance to avoidance 
#                 points weighting, 
#            1 - distance to avoidance points weighing is square-rooted, 
#            2 - distance to avoidance points weighting is squared
# avoidMod  A coefficient to be applied to the avoidance points weighting.
ELE_avoidLocs <- data.frame(
  "x" = c(-155000, -140000, -125000),
  "y" = c(7900000, 7920000, 7900000))
ELE_avoidTransformation <- 2
ELE_avoidModifier <- 4

# additional cycles
ELE_rest_Cycle <- c(0.12, 0, 24, 24)
c0 <- c(0.075, 0, 24* (365/2), 24* 365) # seasonal
ELE_additional_Cycles <- rbind(c0)

# choose start location and options
start <- c(-134900, 7964000)
timesteps <- 1000000 #24*60*31
des_options=10; options=12;

# debug options far from fence
start <- c(-125000, 7920000)
# des_options=2; options=5;
# timesteps = 5


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                         PLOT SIMULATION PARAMETERS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# # params for plots
# move.rast <- terra::rast(procpath('movementRaster.tif'))
# f <- projectMe(fences, 32735)
# b <-1000
# 
# # zoomed out
# plot(move.rast, main="zoomed out")
# plot(f, add=TRUE, col="black")
# points(start[1], start[2], pch=19, col='red')
# 
# # zoomed in
# plot(move.rast, main="zoomed in",
#      xlim=c(start[1]-b, start[1]+b),
#      ylim=c(start[2]-b, start[2]+b))
# plot(f, add=TRUE, col="black")
# points(start[1], start[2], pch=19, col='red', lwd=5)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                         SIMULATE ELEPHANT MOVEMENTS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# roxygen2::roxygenize("~/R_Packages/abmAME")

seed = 10010
set.seed(seed)
# runSim(barriers=barriers, perm=perm,
#        seed=seed, timesteps=timesteps,
#        activities = "all", colorby='inx')
simRes <- abmAME::abm_simulate(
  start = start,
  timesteps = timesteps,
  des_options = des_options,
  options = options,
  shelterLocations = ELE_shelterLocs,
  shelterSize = ELE_shelterSize,
  avoidPoints = ELE_avoidLocs,
  destinationRange = ELE_destinationRange,
  destinationDirection = ELE_destinationDirection,
  destinationTransformation = ELE_destinationTransformation,
  destinationModifier = ELE_destinationModifier,
  avoidTransformation = ELE_avoidTransformation,
  avoidModifier = ELE_avoidModifier,
  k_step = ELE_k_step,
  s_step = ELE_s_step,
  mu_angle = ELE_mu_angle,
  k_angle = ELE_k_angle,
  rescale_step2cell = ELE_rescale,
  behave_Tmat = ELE_behaveMatrix,
  rest_Cycle = ELE_rest_Cycle,
  additional_Cycles = ELE_additional_Cycles,
  landscape_data = ELE_landsdata
)

## SAVING PLOTS
f <- list(projectMe(fences, 32735), projectMe(rivers, 32735))
p <- plotPaths(simRes, f, perm, seed, activity='move', colorby='inx')
filename <- here::here(outdir, "paths", paste0("path_", seed, ".png"))
message('saving plot to: ', filename)
ggsave(filename=filename, plot=p)

## SAVING SIMULATIONS
# filename2 <- here::here(outdir, "simulations", paste0("simulation_", seed, ".rdata"))
# message('saving simulation to: ', filename2)
# save(simRes, file=filename2)
beepr::beep()

# # plotting all
# plot_list <- list(
#   move   = plotPaths(simRes, barriers, perm, seed, activity = 'move',   colorby='inx'),
#   forage = plotPaths(simRes, barriers, perm, seed, activity = 'forage', colorby='inx'),
#   shelter= plotPaths(simRes, barriers, perm, seed, activity = 'shelter',colorby='inx')
# )
# p <- ggpubr::ggarrange(plotlist = plot_list, nrow=1, common.legend=TRUE)
