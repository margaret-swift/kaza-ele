# saveSimulationParameters.R
# Created 16 Oct 2023
# Margaret Swift <mes114@duke.edu>

# Marshall BM and Duthie AB. abmAnimalMovement: An R package for simulating 
#   animal movement using an agent-based model [version 1; peer review: 2 
#   approved with reservations]. F1000Research 2022, 11:1182 
#   (https://doi.org/10.12688/f1000research.124810.1)

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

source(here::here('02_scripts','utilities.R'))
setDataPaths('elephant', verbose=FALSE)
load(procpath('ele.rdata'))

# build abmAME
roxygen2::roxygenize("~/R_Packages/abmAME")
# devtools::install_github("margaret-swift/abmAME")
# library(abmAME)


# ******************************************************************************
#                     Get transition matrix and params from HMM
# ******************************************************************************

# sex-specific parameters
# Real fence encounter rate data from Naidoo et al 2022, at 1km threshold
# for fences, rivers, and roads
perms_f <- c(0, 0.101, 0.153)
perms_m <- c(0.035, 0.145, 0.258)

## Load and/or generate raster and barrier layers
setDataPaths('geographic', verbose=FALSE)
load(procpath('geographic.rdata'))
shelt.rast <- terra::rast(procpath('shelterRaster.tif'))
move.rast <- terra::rast(procpath('movementRaster.tif'))
forage.rast <- terra::rast(procpath('forageRaster.tif'))
barriers <- list(fences, rivers)
rasters <- list(shelt.rast, move.rast, forage.rast)
ELE_landsdata_F <- generateLandscapeData(rasters, barriers, perms_f[1:2])
ELE_landsdata_M <- generateLandscapeData(rasters, barriers, perms_m[1:2])

# ******************************************************************************
#                     Get transition matrix and params from HMM
# ******************************************************************************

createBehaveMatrix <- function(mle) {
  bm <- matrix(nrow=3, ncol=3)
  for (i in 1:3) {
    # get transition matrix data from HMM MLE
    d2 <- exp(mle$beta[,(i*2)-1])
    d3 <- exp(mle$beta[,(i*2)])
    #calculate probs
    p2 <- d2 / ( 1 + d2 + d3 )
    p3 <- d3 / ( 1 + d2 + d3 )
    p1 <- 1 - ( p2 + p3 )
    #save
    bm[i,] <- c( p1, p2, p3 )
  }
  bm
}
load(here(outdir, 'hmmLongF.rdata'))
ELE_behaveMatrix_F <- createBehaveMatrix(hmm.f$mle)
load(here(outdir, 'hmmLongM.rdata'))
ELE_behaveMatrix_M <- createBehaveMatrix(hmm.m$mle)

getSLTA <- function(mle) {
  # step length (CPP CODE ASSUMES MINUTE SCALE)
  k_step <- mle$step['mean',] / 60
  s_step <- mle$step['sd',] / 60
  
  # turning angle
  mu_angle <- mle$angle['mean',]
  k_angle  <- mle$angle['concentration',]
  
  # all together 
  return(list(k_step=k_step, s_step=s_step, mu_angle=mu_angle, k_angle=k_angle))
}

ELE_SLTA_F = getSLTA(hmm.f$mle)
ELE_SLTA_M = getSLTA(hmm.m$mle)


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
ELE_rescale <- ELE_landsdata_F$landscape_meta[5] #should be around 10

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


# ******************************************************************************
## SAVING DATA
filename <- here::here(outdir, "simulations", "simulation_parameters.rdata")
message('saving simulation data to: ', filename)
save(start, timesteps, des_options,options,
     ELE_shelterLocs,ELE_shelterSize, ELE_avoidLocs,
     ELE_destinationRange,ELE_destinationDirection,
     ELE_destinationTransformation,ELE_destinationModifier,
     ELE_avoidTransformation,ELE_avoidModifier,
     ELE_k_step,ELE_s_step,ELE_mu_angle,ELE_k_angle,ELE_rescale,
     ELE_rest_Cycle,ELE_additional_Cycles,
     ELE_behaveMatrix_F,ELE_behaveMatrix_M,
     ELE_landsdata_F,ELE_landsdata_M, 
     ELE_SLTA_F, ELE_SLTA_M,
     file=filename
     )
