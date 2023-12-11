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
               nngeo,
               reshape2,
               ggplot2,
               terra)
source(here('02_scripts','utilities.R'))
setDataPaths('elephant', verbose=FALSE)
load(procpath('ele.rdata'))

setDataPaths('geographic', verbose=FALSE)
load(procpath('geographic.rdata'))
shelt.rast <- terra::rast(procpath('shelterRaster.tif'))
move.rast <- terra::rast(procpath('movementRaster.tif'))
forage.rast <- terra::rast(procpath('forageRaster.tif'))

barriers <- list(fences, rivers)
rasters <- list(shelt.rast, move.rast, forage.rast)

# ******************************************************************************
#                                   FUNCTIONS
# ******************************************************************************
# MATRIX CREATION
.createMat <- function(rast) {
  # transforming raster to matrix with size and extent info
  mat <- as.matrix(rast, wide=TRUE)
  mat[is.nan(mat)] <- -99.9
  mat
}
# CHECK INPUTS
.checkMat <- function(mat) {
  minCheck <- min(mat, na.rm=TRUE) >= -99.9
  maxCheck <- max(mat, na.rm=TRUE) <= 1
  flag = all(minCheck, maxCheck)
  flag
}
# helper functions not to export
.segmentBarriers <- function(b, id=0, ...) {
  data <- st_segments(b, ...)
  NR = nrow(data)
  dig <- floor(log10(NR)) + 1
  data$INX <- id + (10 ^ -dig) * (1:NR)
  return(data)
}
.createMatrix <- function(data, p) {
  coords <- data %>% st_coordinates()
  NR <- nrow(coords)
  # fill in barrier matrix
  b.mat <- matrix(data=0, nrow=NR/2, ncol=6)
  b.mat[,1:2] <- coords[seq(1, NR, by=2), 1:2]
  b.mat[,3:4] <- coords[seq(2, NR, by=2), 1:2]
  b.mat[,5] <- p
  b.mat[,6] <- data$INX
  colnames(b.mat) = c('x', 'y', 'xend', 'yend', 'perm', 'id')
  b.mat
}
# .createLookup <- function(data, rast) {
#   tab <- terra::extract(rast, data)
#   tab$INX <- data$INX[tab$ID]
#   inx.na <- is.na(tab[,2])
#   tab <- tab[!inx.na,]
#   tab <- tab[order(tab[,2]),]
#   tab <- as.matrix(tab)
#   return(tab)
# }
.generateBarriers <- function(barriers, p_list, ...) {
  # inputs are a shapefile of barrier data and a list of accompanying
  # permeabilities for each one. If all are the same for one barrier type,
  # then we only have to run it once for the whole dataset
  LP <- length(p_list)
  if (class(barriers)[1] == 'list') NB <- length(barriers)
  else NB <- nrow(barriers)
  if (LP != NB) {
    message("WARNING! length of permeability vector must be 1 or # of barriers")
    return()
  } else {
    if (LP == 1) p <- p_list
    df = NULL
    for (i in 1:NB) {
      if (LP > 1) p = p_list[i]
      b <- barriers[[i]]
      if (nrow(b) > 0) { # if cropping hasn't removed all data
        segs <- .segmentBarriers(b, i, progress=FALSE)
        mat  <- .createMatrix(segs, p)
        if (is.null(df)) df <- mat
        else df <- rbind(df, mat) 
      }
    }
  }
  return(df)
}

# SETTING UP ENVIRONMENT
rectifyLandscapeFeatures <- function(rasters, barriers, perms=1, crs=32735){
  require(terra)
  require(sf)
  # local functions
  .projectMe <- function(obj) {
    mycrs = crs(obj, describe=T)$code
    if (mycrs != crs) {
      if (class(obj)[1] == "SpatRaster") obj = project(obj, paste0("EPSG:", crs)) 
      if (class(obj)[1] == "sf") obj = st_transform(obj, crs) 
    }
    return(obj)
  }
  
  message("== RASTERS ==") 
  message("reprojecting rasters if necessary...")
  rast.reproj = lapply(rasts, .projectMe)
  
  message("acquiring extent and resolution of rasters...")
  rast <- rast.reproj[[1]]
  ext <- as.list(ext(rast))
  envExt <-c(unlist(ext, use.names=FALSE), res(rast))
  
  message('creating raster matrices for landscape variables...')
  matrices <- lapply(rast.reproj, .createMat) 
  checks <- lapply(matrices, .checkMat) 
  if ( sum(unlist(checks)) != length(checks) ) {
    stop("All the landscape layers should be numeric matricies, with values between -99.9 and 1")
  }
  
  message("== BARRIERS ==")
  message("reprojecting barriers...")
  barr.reproj <- lapply(barriers, .projectMe)
  
  message("cropping barriers to raster extent...")
  .cropMe <- function(x) {
    st_agr(x) = "constant"
    st_crop(x, unlist(ext))
  }
  barr.crop <- sapply(barr.reproj, .cropMe)
  
  message('generating barrier data') 
  barr.df <- .generateBarriers(barr.crop, perms, progress=FALSE)
  
  mydata <- list(
    landscape_data = matrices,
    landscape_meta = envExt,
    barrier_data = barr.df
  )
  return(mydata)
}




# ******************************************************************************
#                                   MAIN
# ******************************************************************************

myData <- rectifyLandscapeFeatures(rasters, barriers, perms, crs=32735)
  
