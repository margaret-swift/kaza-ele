# runSimulation.R
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
source(here::here('02_scripts','utilities.R'))
load(here::here("03_output", "simulations", "simulation_parameters.rdata"))
setDataPaths('geographic', verbose=FALSE)
load(procpath('geographic.rdata'))
move.rast <- terra::rast(procpath('movementRaster.tif'))
shelter.rast <- terra::rast(procpath('shelterRaster.tif'))
f <- projectMe(fences, 32735)

## build abmAME
roxygen2::roxygenize("~/R_Packages/abmAME")
# devtools::install_github("margaret-swift/abmAME")
# library(abmAME)

# set flag for which sex to simulate
sex <- "F" #"M"
if (sex == "F") {
  ELE_behaveMatrix = ELE_behaveMatrix_F
  ELE_landsdata = ELE_landsdata_F
  ELE_SLTA = ELE_SLTA_F
} else { 
  ELE_behaveMatrix = ELE_behaveMatrix_M
  ELE_landsdata = ELE_landsdata_M
  ELE_SLTA = ELE_SLTA_F
}

# choose start location and options
start <- c(-134900, 7964000)
des_options=10; options=12;
ELE_shelterLocs <- data.frame(
  "x" = c(-134700, -134000, -133000, 
          -135500, -135000),
  "y" = c(7963800, 7963600, 7963400,
          7963800, 7964500))

# setting up home range area
hr_area = 5000 %>% units::set_units('ha')
units(hr_area) <- 'm2' # transform units to meters2
home_r = sqrt(as.numeric(hr_area) / pi)
homerange <- data.frame(x=start[1], y=start[2])
home_xy = st_as_sf(homerange, coords=c('x', 'y'), crs=32735)
home_circ <- st_buffer(home_xy, dist=home_r)
s <- spatSample(shelter.rast,
                size = 40,
                method = "weights",
                as.points=TRUE,
                ext=ext(home_circ))
ELE_shelterLocs <- geom(s)[,c('x', 'y')] %>% as.data.frame()
ELE_shelterSize <- 5000

# ELE_shelterLocs <- data.frame(
#   "x" = c(-133000),
#   "y" = c(7963400))
# ELE_shelterSize <- 5000

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                         PLOT SIMULATION PARAMETERS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# # params for plots
# b <-2000
# 
# # zoomed out
# plot(move.rast, main="zoomed out")
# plot(f, add=TRUE, col="black")
# points(start[1], start[2], pch=19, col='red')
# points(ELE_shelterLocs$x, ELE_shelterLocs$y, pch=19, col='blue')
# 
# 
# # zoomed in
# plot(move.rast, main="zoomed in",
#      xlim=c(start[1]-b, start[1]+b),
#      ylim=c(start[2]-b, start[2]+b))
# plot(f, add=TRUE, col="black")
# points(start[1], start[2], pch=19, col='red', lwd=5)
# points(ELE_shelterLocs$x, ELE_shelterLocs$y, pch=19, col='blue')


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                         SIMULATE ELEPHANT MOVEMENTS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

seed = 888902

timesteps <- 5000000 #24*60*31
set.seed(seed)
message('seed: ', seed)
# roxygen2::roxygenize("~/R_Packages/abmAME")

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
  home_range_area = hr_area,
  home_range_units = "ha",
  avoidPoints = ELE_avoidLocs,
  destinationRange = ELE_destinationRange,
  destinationDirection = ELE_destinationDirection,
  destinationTransformation = ELE_destinationTransformation,
  destinationModifier = ELE_destinationModifier,
  avoidTransformation = ELE_avoidTransformation,
  avoidModifier = ELE_avoidModifier,
  k_step = ELE_SLTA$k_step,
  s_step = ELE_SLTA$s_step,
  mu_angle = ELE_SLTA$mu_angle,
  k_angle = ELE_SLTA$k_angle,
  rescale_step2cell = ELE_rescale,
  behave_Tmat = ELE_behaveMatrix,
  rest_Cycle = ELE_rest_Cycle,
  additional_Cycles = ELE_additional_Cycles,
  landscape_data = ELE_landsdata
)


## SAVING PLOTS
p <- plotPaths(simRes=simRes, barriers=f, 
               perm=NULL, seed=seed, activity='move', colorby='inx', b=1000)
filename <- here::here(outdir, "paths", paste0("path_", seed, ".png"))
message('saving plot to: ', filename)
ggsave(filename=filename, plot=p, width=10, height=10)

## SAVING SIMULATIONS
# filename <- here::here(outdir, "simulations", paste0("simulation_", seed, ".rdata"))
# message('saving simulation to: ', filename)
# save(simRes, file=filename)
# beepr::beep()

# # plotting all
# plot_list <- list(
#   move   = plotPaths(simRes, barriers, perm, seed, activity = 'move',   colorby='inx'),
#   forage = plotPaths(simRes, barriers, perm, seed, activity = 'forage', colorby='inx'),
#   shelter= plotPaths(simRes, barriers, perm, seed, activity = 'shelter',colorby='inx')
# )
# p <- ggpubr::ggarrange(plotlist = plot_list, nrow=1, common.legend=TRUE)
