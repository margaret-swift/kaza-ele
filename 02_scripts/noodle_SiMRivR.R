# noodleSiMRiv.R
# Created 15 Sept 2023
# Margaret Swift <mes114@duke.edu>

# Noodle file for messing around with SiMRiv:
# 
# Quaglietta, L. & Porto, M. SiMRiv: an R package for mechanistic simulation 
#     of individual, spatially-explicit multistate movements in rivers, 
#     heterogeneous and homogeneous spaces incorporating landscape bias. 
#     Mov Ecol 7, 11 (2019).
# GitHub page: https://github.com/miguel-porto/SiMRiv
# Vignette: https://cran.r-project.org/web/packages/SiMRiv/vignettes/SiMRiv.pdf

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('./02_scripts/noodle_SiMRivR.R')
source('./02_scripts/utilities.R')
p_load("SiMRiv", #v1.0.6
       "raster") 


# ******************************************************************************
#                           MESSING AROUND WITH SIMRIV
# ******************************************************************************

# define a species with a single-state movement type
# characterized by a random walk
SL = 15
rand.walker <- species(state.CRW(0.98) + SL)
sim.rw <- simulate(rand.walker, 10000)
plot(sim.rw, type = "l", asp = 1, main = "Random walk")


# Levy walks
tm <- transitionMatrix(0.005, 0.01, #state 1>2 and 1>3
                       0.05,  0.02, #state 2>1 and 2>3
                       0.05,  0.005)#state 3>1 and 3>2
levy.walker <- species(
  state.RW() + state.CRW(0.98) + state.Resting(), 
  trans = tm) + SL
sim.lw <- simulate(levy.walker, 1e5)
plot(sim.lw, type='l', asp=1, main='levy walker')



# ******************************************************************************
#                           MESSING AROUND WITH RASTERS
# ******************************************************************************

# import some data
setDataPaths('geographic')
load(here(procpath, 'fences.Rdata')) #two fences: 'fence' and 'cbpp'
load(here(procpath, 'geographicData.RData'))

# set resistances for different landscape types
# 0 = no resistance
# 1 = total resistance
lands.meta$resistance = c(
  0, # bare area
  0, # bare floodplain area
  1, # built-up
  0.7, # closed bushland
  0.9, # closed forest
  0.7, # closed herbaceous wetland
  0.9, # closed woodland
  0.2, # cropland
  0.3, # open bushland/shrubs
  0.1, # open herbaceous vegetation
  0.1, # open herbaceous wetland
  0.1, # open herbaceous floodplain
  0.2, # open woodland/bushland
  0.2, # sparse forest/woodland
  0.1, # sparse herbaceous wetland
  0.2, # sparse/open bushland/shrubs
  1, # water bodies permanent
  0.7  # water bodies seasonal
)

x1 <- reclassify(lands.raster.1, lands.meta[,c('category', 'resistance')])
x2 <- reclassify(lands.raster.2, lands.meta[,c('category', 'resistance')])
crs(x1) <- crs(x2) <- crs(fence)
e1 <- extent(x1)
e2 <- extent(x2)
e <- merge(e1, e2)

x11 <- brick(x1)
x22 <- brick(x2)
extent(x11) <- extent(x22) <- e
s <- mosaic(x11, x22, tolerance=0.5, fun=sum)


land.files <- list.files(rawpath, full.names=TRUE)

# read LC categories metadata
meta.file <-  land.files[grepl('csv', land.files)]
lands.meta <- read.csv(meta.file)

# read shapefile
shp.file <- land.files[grepl('tif', land.files)]
lc.r <- raster(shp.file)

# crop LC to a reprojected khaudum map, then reproject smaller raster to khau crs

bb1 <- bb2 <- st_bbox(collar.df) 
bb1 [['ymin']] <- bb2 [['ymax']] <- -18.59715

bb1 <- bb1 %>% 
  nominatimlite::bbox_to_poly() %>% 
  st_transform(crs(lc.r)) %>% 
  st_as_sf()

bb2 <- bb2 %>% 
  nominatimlite::bbox_to_poly() %>% 
  st_transform(crs(lc.r)) %>% 
  st_as_sf()

t0 <- proc.time()[['elapsed']]
lands.raster.1 <- crop(lc.r, bb1) %>%
  projectRaster(crs=crs(khau), method="ngb") #### THIS TAKES FOREVER FYI ! ####
lands.raster.2 <- crop(lc.r, bb2)%>%
  projectRaster(crs=crs(khau), method="ngb")
lands.raster <- stack(lands.raster.1, lands.raster.2)
t <- proc.time()[['elapsed']] - t0
message(round(t/60, 2), " minutes elapsed")

# now for the fences
# helpful: https://michaelbcalles.netlify.app/post/rasterize-lines-in-r/
temp <- raster(e, resolution = res(x1), crs = crs(x1))
y <- resistanceFromShape(st_union(fence, cbpp), # function from simriv
                         baseRaster=temp,
                         binary=TRUE,
                         field=1,
                         background=0,
                         buffer=50,
                         extend=FALSE)
y <- reclassify(y, cbind(0, 1))

cols <- hcl.colors(12, "YlOrRd", rev = TRUE)
breaks <- seq(0, 1, length.out=13)
image(y, col=cols, breaks=breaks)
image(x1, add=TRUE, col=cols, breaks=breaks)
image(x2, add=TRUE, col=cols, breaks=breaks)
image(y, add=TRUE, col=cols, breaks=breaks)

x = merge(x1, x2, tolerance=1)




