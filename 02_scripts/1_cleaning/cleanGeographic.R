# cleanGeographic.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanGeographic.R')
source(here::here('02_scripts', 'utilities.R'))

# ******************************************************************************
# load boundaries, linear, and water features and save in group files
# these can then be loaded anywhere using 'quickload()':
# quickload()

setDataPaths('boundaries')
kaza = st_read(rawpath('kaza_boundary'))
kaza_aoi = st_read(rawpath('kaza_aoi'))
khau = st_read(rawpath('khaudum_boundary'))
khau_tall_aoi = st_read(rawpath('khau_tall_aoi'))
namib_ele_aoi = st_read(rawpath('namib_ele_aoi'))

ggplot() + 
  geom_sf(data = kaza, color='green', fill=NA) +
  geom_sf(data = khau, color='purple', fill=NA) +
  geom_sf(data = kaza_aoi, color='orange', fill=NA) +
  geom_sf(data = khau_tall_aoi, color='blue', fill=NA) +
  geom_sf(data = namib_ele_aoi, color='black', fill=NA)

save(kaza, kaza_aoi, khau, khau_tall_aoi, namib_ele_aoi, 
     file=procpath('boundaries.rdata'))


# ******************************************************************************
# load linear features and save to one rdata file.
setDataPaths('linear_features')
fences = st_read(rawpath('merged_fences_old'))
rivers = st_read(rawpath('rivers_main_KAZA'))
roads = st_read(rawpath('primary_trunk_rds'))

ggplot() + 
  geom_sf(data = fences, color='green', fill=NA) +
  geom_sf(data = rivers, color='purple', fill=NA) +
  geom_sf(data = roads, color='orange', fill=NA)

save(fences, rivers, roads, file=procpath('linear_features.rdata'))

# ******************************************************************************
# load waters and save to one rdata file.
setDataPaths('waters')
waters_art = st_read(rawpath('khaudum_artificial_waters'))
waters_nat = st_read(rawpath('waterhole_stats'))
save(waters_art, waters_nat, file=procpath('waters.rdata'))





# ******************************************************************************
#                         Arthur's fence status data
# ******************************************************************************

# setDataPaths('linear_features')
# library('readxl')
# pth <- "fence_damage_arthur"
# files <- list.files(rawpath(pth))
# 
# ZF_locs <- read_excel(rawpath(pth, "Zambezi.Notes.fencecondition.Dec2020.xlsx"), col_names=TRUE)
# WBF_locs <- read_excel(rawpath(pth, "WBF.5km.coords.xlsx"), col_names =FALSE)
# names(WBF_locs)

# ******************************************************************************
#                        Set up raster resistance layer
# ******************************************************************************


setDataPaths('landcover')
lands.meta = read.csv(rawpath('landcover_metadata.csv'))
lr <- terra::rast(rawpath('eoss4wwf_kaza_tfca_landcover_2020.tif'))

# kinda the only way I could find to mosaic rasters with different extents...
# https://stackoverflow.com/questions/67169266/error-in-mosaic-of-rasters-from-different-extent-using-terra-package-in-r
# f1 <- here('03_output', 'tmp', 'test1.tif')
# f2 <- here('03_output', 'tmp', 'test2.tif')
# r1 <- writeRaster(lands.raster.1, f1, overwrite=TRUE)
# r2 <- writeRaster(lands.raster.2, f2, overwrite=TRUE)
# lands.vrt <- terra::vrt(c(f1, f2), "03_output/tmp/test.vrt", overwrite=TRUE) 

# define generic landcover type bins
lands.meta$cover = c(
  0, # bare area
  0, # bare floodplain area
  0, # built-up
  2, # closed bushland
  2, # closed forest
  2, # closed herbaceous wetland
  2, # closed woodland
  1, # cropland
  3, # open bushland/shrubs
  3, # open herbaceous vegetation
  3, # open herbaceous wetland
  3, # open herbaceous floodplain
  3, # open woodland/bushland
  4, # sparse forest/woodland
  4, # sparse herbaceous wetland
  4, # sparse/open bushland/shrubs
  0,   # water bodies permanent
  0    # water bodies seasonal
)
cover_class = data.frame(cover=0:4, 
                         cover_class=c('bare/water', 'cropland', 'closed', 
                                       'open', 'sparse'))
lands.meta = left_join(lands.meta, cover_class, by="cover")

# define generic vegetation type bins
lands.meta$veg = c(
  0, # bare area
  0, # bare floodplain area
  0, # built-up
  2, # closed bushland
  3, # closed forest
  4, # closed herbaceous wetland
  3, # closed woodland
  1, # cropland
  2, # open bushland/shrubs
  4, # open herbaceous vegetation
  4, # open herbaceous wetland
  4, # open herbaceous floodplain
  3, # open woodland/bushland
  3, # sparse forest/woodland
  4, # sparse herbaceous wetland
  2, # sparse/open bushland/shrubs
  0,   # water bodies permanent
  0    # water bodies seasonal
)
veg_class = data.frame(veg=0:4, 
                       veg_class=c('nonveg', 'cropland', 'bushland', 
                                   'forest/woodland', 'herbaceous/wet'))
lands.meta = left_join(lands.meta, veg_class, by="veg")

# ## RECLASSIFY RASTERS FOR ACTIVITY STATES
# st_transform(khau, crs="EPSG:32734")
# bb <- terra::ext(c(400000, 600000, 7200000, 8300000))
# lr.crop <- terra::crop(lr, bb)
# plot(lr.crop)
# reclassWrite <- function(type) {
#   reclass = terra::classify(lr.crop, lands.meta[,c('category', type)])
#   reproj = terra::project(reclass, 'EPSG:32735')
#   terra::writeRaster(reproj,
#               filename=procpath(paste0(type, 'Raster.tif')),
#               overwrite=TRUE)
# }

write.csv(lands.meta, file=rawpath('kaza_landcover', 'landcover_metadata.csv'))
