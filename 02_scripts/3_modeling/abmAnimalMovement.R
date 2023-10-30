# abmAnimalMovement.R
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
load(here(outdir, 'hmmLongF.rdata'))
load(here(outdir, 'hmmLongM.rdata'))

setDataPaths('geographic')
load(procpath('geographic.rdata'))
# location for updated abmAnimalMovementMES
pkg.path <- here("02_scripts", "9_packages", "abmAnimalMovementMES")
# devtools::load_all(path=pkg.path)


# ******************************************************************************
#                         Get transition matrix from HMM
# ******************************************************************************
# 
# # get estimates for females
# hmm <- hmm.f
# 
# vals <- 1:3
# probs.df <- data.frame(matrix(nrow=1, ncol=9))
# names(probs.df) <- paste(rep(vals, each=3), vals, sep=" -> ")
# mat <- hmm$mle$beta
# for (i in vals) {
#   types <- paste(i, vals[vals != i], sep=" -> ")
# 
#   d2 <- exp(mat[,types[1]])
#   d3 <- exp(mat[,types[2]])
#   sum <- d2 + d3
# 
#   p2 <- d2 / ( 1 + sum )
#   p3 <- d3 / ( 1 + sum )
#   p1 <- 1 - ( p2 + p3 )
# 
#   probs <- data.frame( p1, p2, p3 )
#   names(probs) <- c(paste(i, i, sep=" -> "), types)
#   probs.df[,names(probs)] <- probs
# }
# probs.df <- unlist(probs.df)
probs.df <- c(0.68719185, 0.2962346, 0.01657354,
  0.13275791, 0.7859501, 0.08129195,
  0.01519633, 0.2423039, 0.74249974)
tm <- matrix(probs.df, nrow=3, byrow=TRUE)


# ******************************************************************************
#                       Trying out with the badger example
# ******************************************************************************

row <- 200; col <- 200;
landscapeLayersList <- list(
  "shelter" = matrix(runif(row*col, 0, 1), nrow = row, ncol = col),
  "forage" = matrix(runif(row*col, 0, 1), nrow = row, ncol = col),
  "movement" = matrix(runif(row*col, 0, 1), nrow = row, ncol = col))
BADGER_shelter <- landscapeLayersList$shelter
BADGER_forage <- landscapeLayersList$forage
BADGER_move <- landscapeLayersList$movement

# # give them only a few good foraging spots
# rast <- BADGER_forage * 0 + runif(n=row*col, 0, 0.25)
# rast[80:85,80:85] <- 1
# BADGER_forage <- rast






# shelters
BADGER_shelterLocs <- data.frame(
  "x" = c(90, 95, 110),
  "y" = c(90, 75, 90))
BADGER_shelterSize <- 5

# b0 <- c(0.97, 0.01, 0.001) # rest
# b1 <- c(0.0002, 0.95, 0.0008) # explore/move
# b2 <- c(0.001, 0.00001, 0.99) # forage
# BADGER_behaveMatrix <- rbind(b0, b1, b2)
BADGER_behaveMatrix <- tm

# step params
BADGER_k_step <- c(0.3*60, 1.25*60, 0.25*60)
BADGER_s_step <- c(0.8, 0.25, 0.5)
BADGER_mu_angle <- c(0, 0, 0)
BADGER_k_angle <- c(0.6, 0.99, 0.6)

# destination params
BADGER_destinationRange <- c(3, 120)
BADGER_destinationDirection <- c(0, 0.01)
BADGER_destinationTransformation <- 2
BADGER_destinationModifier <- 2

BADGER_rescale <- 50

# avoiding
# avoidTrans 0 - no transformation applied to the distance to avoidance 
#                 points weighting, 
#            1 - distance to avoidance points weighing is square-rooted, 
#            2 - distance to avoidance points weighting is squared
# avoidMod  A coefficient to be applied to the avoidance points weighting.
nrep=1e3
BADGER_avoidLocs <- data.frame(
  "x" = 5,
  "y" = 5)
BADGER_avoidTransformation <- 2
BADGER_avoidModifier <- 4

# additional cycles
BADGER_rest_Cycle <- c(0.12, 0, 24, 24)
c0 <- c(0.075, 0, 24* (365/2), 24* 365) # seasonal
BADGER_additional_Cycles <- rbind(c0)

# startLocation <- sample(90:110, 2, replace = TRUE)
start <- c(95, 95)
timesteps <- 5000#24*60*31
des_options=10
options=12

# Declaring the fence 
fence_points = matrix(data=c(100, 80, 100, 100,
                             100, 80, 115, 80),
                      nrow=2,
                      byrow=TRUE)
fence_display = as.data.frame(fence_points)
names(fence_display) = c('x', 'y', 'xend', 'yend')

# fences.df <- fences[102:103,] %>% 
#   st_cast('MULTIPOINT') 
# ggplot() + 
#   geom_sf(data=fences[102:103,], 
#           color='gray', linewidth=2) +
#   geom_sf(data=fences.df, 
#           mapping=aes(color=Name))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
createMoveRast <- function(simRes) {
  locs <- simRes$locations
  df <- data.frame( x = rep(1:200, each=200),
                    y = rep(1:200, times=200),
                    shelter=NA, forage=NA, move=NA) %>%
    dplyr::filter(x > min(locs$x), x < max(locs$x),
                  y > min(locs$y), y < max(locs$y))
  df
}
myPlot <- function(simRes, seed=1001, inx=NULL) {
  # get rasters 
  move_rast <- createMoveRast(simRes)
  for(i in 1:nrow(move_rast)) {
    y <- move_rast$y[i]; x <- move_rast$x[i];
    move_rast$shelter[i] <- BADGER_shelter[y,x]
    move_rast$forage[i] <- BADGER_forage[y,x]
    move_rast$move[i] <- BADGER_move[y,x]
  }
  
  # get paths of walking
  path <- simRes$locations
  if (!is.null(inx)) path <- path[inx,]
  path$INX <- 1:nrow(path)
  xlim = c(min(path$x), max(path$x))
  ylim = c(min(path$y), max(path$y))
  title = paste0("Ele movement, seed = ", seed)
  p <- ggplot() +
    # geom_tile(data = move_rast,
              # mapping=aes(x=x, y=y, fill=forage),
              # alpha=0.4) +
    geom_segment(data=fence_display, aes(x=x, y=y, xend=xend, yend=yend), linewidth=2) +
    geom_point(data=BADGER_shelterLocs, 
               mapping=aes(x=x, y=y), color="#dbafed", alpha=0.3,
               size=BADGER_shelterSize*19.5) +
    geom_point(data=simRes$locations[1,], 
               mapping=aes(x=x, y=y), color="red", alpha=0.3, 
               size=20) +
    geom_path(data=path,
              mapping=aes(x=x, y=y, color=INX),
              linewidth=1) +
    geom_point(data=simRes$locations[1,],
               mapping=aes(x=x, y=y),
               color='black', size=5) +
    geom_point(data=simRes$locations[1,],
               mapping=aes(x=x, y=y),
               color='white', size=3) +
    geom_point(data=BADGER_shelterLocs, 
               mapping=aes(x=x, y=y), color="purple", size=4) +
    scale_color_distiller(palette='Spectral', direction=1) +
    scale_fill_distiller(palette="Greens", direction=1) +
    coord_sf(xlim=xlim, ylim=ylim) +
    ggtitle(title)
  p
}

runSim <- function(seed=1001) {
  # Good seeds: 292, 
  set.seed(seed)
  simRes <- abmAnimalMovementMES::abm_simulate(
    start = start, timesteps = timesteps, des_options = des_options,options = options,
    shelterLocations = BADGER_shelterLocs,shelterSize = BADGER_shelterSize,
    avoidPoints = BADGER_avoidLocs,destinationRange = BADGER_destinationRange,
    destinationDirection = BADGER_destinationDirection,
    destinationTransformation = BADGER_destinationTransformation,
    destinationModifier = BADGER_destinationModifier, avoidTransformation = BADGER_avoidTransformation,
    avoidModifier = BADGER_avoidModifier,k_step = BADGER_k_step,s_step = BADGER_s_step,
    mu_angle = BADGER_mu_angle,k_angle = BADGER_k_angle, rescale_step2cell = BADGER_rescale,
    behave_Tmat = BADGER_behaveMatrix,rest_Cycle = BADGER_rest_Cycle,
    additional_Cycles = BADGER_additional_Cycles,shelteringMatrix = BADGER_shelter,
    foragingMatrix = BADGER_forage,movementMatrix = BADGER_move,
    fence=fence_display
  )
  message('seed ', seed)
  p <- myPlot(simRes, seed)
  print(p)
}

# devtools::load_all(path=pkg.path)'

randseeds <- floor(runif(1e5, 0, 1) * 1e5)
i = 0

# run this as a chunk
# # good seeds: 85984; 61862; 78735; 50575
i = i+1
seed <- randseeds[i]
runSim(seed)



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++












# 
# 
# 
# 
# 
# # my visualizations
# move_rast <- data.frame(x = rep(1:200, each=200),
#                         y = rep(1:200, times=200),
#                         shelter=NA, forage=NA, move=NA)
# filterToData <- function(df) {
#   locs <- simRes$locations
#   df <- df %>%
#     dplyr::filter(x > min(locs$x),
#            x < max(locs$x),
#            y > min(locs$y),
#            y < max(locs$y),
#     )
#   df
# }
# move_rast <- filterToData(move_rast)
# for(i in 1:nrow(move_rast)) {
#   y <- move_rast$y[i]; x <- move_rast$x[i];
#   move_rast$shelter[i] <- BADGER_shelter[y,x]
#   move_rast$forage[i] <- BADGER_forage[y,x]
#   move_rast$move[i] <- BADGER_move[y,x]
# }
# 
# 
# myPlot <- function(base = "move", inx=1:nrow(simRes$locations)) {
#   # dat <- move_rast %>%
#   #   dplyr::select_at(c('x', 'y', base))
#   # names(dat) <- c('x', 'y', 'base')
#   # dat$base <- factor(dat$base)
#   path <- simRes$locations[inx,]
#   path$INX <- 1:nrow(path)
#   xlim = c(min(path$x), max(path$x))
#   ylim = c(min(path$y), max(path$y))
#   title = paste0("Ele movement, seed = ", seed)
#   p <- ggplot() +
#     geom_rect(aes(xmin=100, xmax=105, ymin=-Inf, ymax=Inf), fill='darkgray') +
#     geom_segment(data=fence_display, aes(x=x, y=y, xend=xend, yend=yend), linewidth=2) +
#     geom_path(data=path,
#               mapping=aes(x=x, y=y, color=INX),
#               linewidth=1) +
#     geom_point(data=simRes$locations[1,],
#                mapping=aes(x=x, y=y),
#                color='black', size=5) +
#     geom_point(data=simRes$locations[1,],
#                mapping=aes(x=x, y=y),
#                color='white', size=3) +
#     scale_color_distiller(palette='Spectral', direction=1) +
#     coord_sf(xlim=xlim, ylim=ylim) +
#     ggtitle(title)
#   p
# }
# 
# 
# myPlot()

# p2 <- myPlot('forage')
# p3 <- myPlot('shelter')
# p1 + p2 + p3 + plot_layout(guides='collect', nrow=1)



simRes <- abmAnimalMovementMES::abm_simulate(
  # a data frame with x and y coordinates
  start = start,
  # an integer describing the length of the simulation
  timesteps = timesteps,
  # an integer describing the number of foraging destination options an animal
  # is offered
  des_options = des_options,
  # an integer describing the number of movement options an animal is offered
  options = options,
  # a data frame providing x and y coordinates of the shelter locations
  shelterLocations = BADGER_shelterLocs,
  # a value describing the radius around shelter sites that movement step
  # lengths are reduced
  shelterSize = BADGER_shelterSize,
  # a data frame providing x and y coordinates of the avoidance locations
  avoidPoints = BADGER_avoidLocs,
  # a numeric vector of length two that contains the shape and scale values
  # that describe the Gamma distribution for potential foraging destinations
  destinationRange = BADGER_destinationRange,
  # a numeric vector of length two that contains the mean and concentration values
  # that describe the Von Mises distribution for potential foraging destinations
  destinationDirection = BADGER_destinationDirection,
  # a value to chose the type of transformation applied to the animal’s
  # attraction to a chosen destination
  destinationTransformation = BADGER_destinationTransformation,
  # a value modifying the animal’s attraction to a chosen destination
  destinationModifier = BADGER_destinationModifier,
  # a value to chose the type of transformation applied to the animal’s
  # avoidance to a avoidance locations
  avoidTransformation = BADGER_avoidTransformation,
  # a value modifying the animal’s avoidance of avoidance locations
  avoidModifier = BADGER_avoidModifier,
  # a vector of three numbers describing the three behavioural state’s Gamma
  # distributions’ shape parameter for step lengths
  k_step = BADGER_k_step,
  # a vector of three numbers describing the three behavioural state’s Gamma
  # distributions’ scale parameter for step lengths
  s_step = BADGER_s_step,
  # a vector of three numbers describing the three behavioural state’s Von Mises
  # distributions’ mean parameter for turn angles
  mu_angle = BADGER_mu_angle,
  # a vector of three numbers describing the three behavioural state’s Von Mises
  # distributions’ concentration parameter for turn angles
  k_angle = BADGER_k_angle,
  # a numeric value to specify the size of the environmental matrices cells
  rescale_step2cell = BADGER_rescale,
  # a 3x3 numeric matrix describing the transition probabilities between the
  # three behavioural states
  behave_Tmat = BADGER_behaveMatrix,
  # A vector length 4 for amplitude, midline, offset and frequency to define the
  # resting/active cycle
  rest_Cycle = BADGER_rest_Cycle,
  # A data.frame 4 columns wide for amplitude, midline, offset and frequency to
  # define any additional activity cycles, where each row is another cycle
  additional_Cycles = BADGER_additional_Cycles,
  # Three arguments for the three matrices describing the landscape the
  # simulated animal occupies
  shelteringMatrix = BADGER_shelter,
  foragingMatrix   = BADGER_forage,
  movementMatrix   = BADGER_move
)

