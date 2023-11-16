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
load(here(outdir, 'hmmLongF.rdata'))
load(here(outdir, 'hmmLongM.rdata'))

setDataPaths('geographic')
load(procpath('geographic.rdata'))

# location for updated abmFences
pkg.path <- "C:/Users/mes473/OneDrive - Cornell University/Documents/R_Packages/abmFences"
roxygen2::roxygenize(package.dir=pkg.path)
# devtools::install_github("margaret-swift/abmFences")
# library(abmFences)

# ******************************************************************************
#                         Get transition matrix from HMM
# ******************************************************************************

# get estimates for females
hmm <- hmm.f


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
probs.df <- unlist(probs.df)
tm <- matrix(probs.df, nrow=3, byrow=TRUE)


# ******************************************************************************
#                             Badger example parameters
# ******************************************************************************

row <- 200; col <- 200;
landscapeLayersList <- list(
  "shelter" = matrix(runif(row*col, 0, 1), nrow = row, ncol = col),
  "forage" = matrix(runif(row*col, 0, 1), nrow = row, ncol = col),
  "movement" = matrix(runif(row*col, 0, 1), nrow = row, ncol = col))
ELE_shelter <- landscapeLayersList$shelter
ELE_forage <- landscapeLayersList$forage
ELE_move <- landscapeLayersList$movement

# shelters
ELE_shelterLocs <- data.frame(
  "x" = c(90, 95, 110),
  "y" = c(90, 75, 90))
ELE_shelterSize <- 5

#behavior matrix
ELE_behaveMatrix <- tm

# step params
ELE_k_step <- c(0.3*60, 1.25*60, 0.25*60)
ELE_s_step <- c(0.8, 0.25, 0.5)
ELE_mu_angle <- c(0, 0, 0)
ELE_k_angle <- c(0.6, 0.99, 0.6)

# destination params
ELE_destinationRange <- c(3, 120)
ELE_destinationDirection <- c(0, 0.01)
ELE_destinationTransformation <- 2
ELE_destinationModifier <- 2

ELE_rescale <- 50

# avoiding
# avoidTrans 0 - no transformation applied to the distance to avoidance 
#                 points weighting, 
#            1 - distance to avoidance points weighing is square-rooted, 
#            2 - distance to avoidance points weighting is squared
# avoidMod  A coefficient to be applied to the avoidance points weighting.
nrep=1e3
ELE_avoidLocs <- data.frame(
  "x" = 5,
  "y" = 5)
ELE_avoidTransformation <- 2
ELE_avoidModifier <- 4

# additional cycles
ELE_rest_Cycle <- c(0.12, 0, 24, 24)
c0 <- c(0.075, 0, 24* (365/2), 24* 365) # seasonal
ELE_additional_Cycles <- rbind(c0)

# startLocation <- sample(90:110, 2, replace = TRUE)
start <- c(95, 95)
timesteps <- 5000#24*60*31
des_options=10
options=12

# # Fake fence data
# fence_points = matrix(data=c(100, 80, 100, 100,
#                              100, 50, 100, 77,
#                              100, 80, 115, 80),
#                       ncol=4,
#                       byrow=TRUE)
# fence_display = as.data.frame(fence_points) %>% 
#   mutate(perm=c(0,0,0.25), id=0)
# names(fence_display) = c('x', 'y', 'xend', 'yend', 'perm', 'id')

# Real fence encounter rate data from Naidoo et al 2022, at 1km encounter threshold
fence_data_f = generateBarriers(list(fences, rivers, roads),
                                   c(0, 0.101, 0.153))
fence_data_m = generateBarriers(list(fences, rivers, roads),
                                   c(0.035, 0.145, 0.258))


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                         SIMULATE ELEPHANT MOVEMENTS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# give them only a few good foraging spots
rast <- ELE_forage * 0 + rbeta(n=row*col, 0.25, 0.25)
rast[85:95, 85:95] <- 1
rast[83:86, 105:110] <- 1
ELE_forage <- rast

# create a movement corridor
rast <- ELE_forage * 0 + runif(n=row*col, 0, 0.1)
rast[0:200, 92:95] <- 1
ELE_move <- rast

timesteps=3e4

## Run these just once at the start of simulations (or when you've updated ++)
# devtools::load_all(path=pkg.path)
set.seed(1001)
randseeds <- floor(runif(1e5, 0, 1) * 1e5)
i = 0

## Run this as a chunk to generate random maps of elephant movements
i = i+1
seed <- randseeds[i]
runSim(seed, colorby='inx')

## Particular seeds that show off the movements well
# - 61862; 78735; 50575; 85984 show fence behaviors really well 
# - 46536; 33712; has good action on both sides of the fence

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




