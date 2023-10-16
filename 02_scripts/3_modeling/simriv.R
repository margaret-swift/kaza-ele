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
       "raster", 'terra') 
load(here(outdir, 'hmmLongF.rdata'))
load(here(outdir, 'hmmLongM.rdata'))
# import fence data
setDataPaths('geographic')
load(here(procpath, 'geographic.rdata')) #two fences: 'fence' and 'cbpp'

v <- vect(fences)
r <- rast(v, nrows=200, ncols=500)
x <- rasterizeGeom(v, r, "crosses")    

m <- c(0, 0.5, 0, 0.5, 4, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- raster( classify(x, rclmat) )
# e <- ext(20.5, 21.5, -19, -18)
# rc <- raster(crop(rc, e))
plot(rc)
lines(v)

# ******************************************************************************
#                           MESSING AROUND WITH SIMRIV
# ******************************************************************************

# Levy walks
hmm <- hmm.f

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
probs.df <- unlist(probs.df)
tm <- matrix(probs.df, nrow=3, byrow=TRUE)
  
# create levy walker
SL = 0.001
levy.walker <- species(
  state.RW() + state.CRW(0.98) + state.Resting(), 
  trans = tm, 
  name="Females") + SL

# initialize point
# init = xyFromCell(rc, sample(which(values(rc) == 0), 1))
init = c(21.8, -18.0)

# start simulation
# seed = seed + 1
seed = 1029
set.seed(seed)
sim.lw <- simulate(individuals = levy.walker, 
                   time = 1e5, 
                   resist=rc,
                   start.resistance=0.5,
                   coords=init)
plot(rc,
     xlim=c(20.5, 22.5),
     ylim=c(-19,-17),
     # ylim = range(sim.lw[, 2]),
     # xlim = range(sim.lw[, 1]),
     main=seed)
# plot(v, col='red')
lines(sim.lw, type='l')
lines(v, lwd=2, col='red')
points(init, cex=1, col='red', pch=19)
