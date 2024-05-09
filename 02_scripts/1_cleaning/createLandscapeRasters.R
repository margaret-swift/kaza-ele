# createLandscapeRasters.R
# Created 04 March 2024
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanGeographic.R')
source(here::here('02_scripts', 'utilities.R'))
setDataPaths('geographic')
load(procpath('geographic.rdata'))
lands.meta = read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))
tfile <- rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif')
lr <- terra::rast(tfile)


library(gdalraster)

# LANDFIRE Existing Vegetation Type (EVT)
evt_file <- system.file(tfile, package="gdalraster")

# make a copy to modify
f <- paste0(tempdir(), "/", "storml_evt_tmp.tif")
file.copy(evt_file,  f)
#> [1] TRUE
ds <- new(GDALRaster, f, read_only=FALSE)
ds$getDefaultRAT(band=1)
pacman::p_load(gdalraster)
tbl <- buildRAT(lr,
                table_type = "thematic",
                na_value = -9999,
                join_df = lands.meta)

# define movement ease 
# 0 = no movement
# 1 = free movement
lands.meta$movement = c(
  1,   # bare area
  1,   # bare floodplain area
  0.1, # built-up
  0.3, # closed bushland
  0.1, # closed forest
  0.3, # closed herbaceous wetland
  0.1, # closed woodland
  0.8, # cropland
  0.7, # open bushland/shrubs
  0.9, # open herbaceous vegetation
  0.9, # open herbaceous wetland
  0.9, # open herbaceous floodplain
  0.8, # open woodland/bushland
  0.8, # sparse forest/woodland
  0.9, # sparse herbaceous wetland
  0.8, # sparse/open bushland/shrubs
  0,   # water bodies permanent
  0.7  # water bodies seasonal
)

# define forage preferences
# 0 = no forage value
# 1 = high forage value
lands.meta$forage = c(
  0, # bare area
  0, # bare floodplain area
  0, # built-up
  0.7, # closed bushland
  0.9, # closed forest
  0.7, # closed herbaceous wetland
  0.9, # closed woodland
  0.9, # cropland
  0.9, # open bushland/shrubs
  0.9, # open herbaceous vegetation
  0.5, # open herbaceous wetland
  0.5, # open herbaceous floodplain
  0.9, # open woodland/bushland
  0.9, # sparse forest/woodland
  0.9, # sparse herbaceous wetland
  0.9, # sparse/open bushland/shrubs
  0,   # water bodies permanent
  0    # water bodies seasonal
)

# define shelter preferences
# 0 = no shelter value
# 1 = high shelter value
lands.meta$shelter = c(
  0, # bare area
  0, # bare floodplain area
  0, # built-up
  0.7, # closed bushland
  0.9, # closed forest
  0.7, # closed herbaceous wetland
  0.9, # closed woodland
  0.9, # cropland
  0.9, # open bushland/shrubs
  0.9, # open herbaceous vegetation
  0.5, # open herbaceous wetland
  0.5, # open herbaceous floodplain
  0.9, # open woodland/bushland
  0.9, # sparse forest/woodland
  0.9, # sparse herbaceous wetland
  0.9, # sparse/open bushland/shrubs
  1,   # water bodies permanent
  1    # water bodies seasonal
)

# COVER CLASS
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

# VEGETATION CLASS
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

## RECLASSIFY RASTERS FOR ACTIVITY STATES
st_transform(khau, crs="EPSG:32734")
bb <- terra::ext(c(400000, 600000, 7200000, 8300000))
lr.crop <- terra::crop(lr, bb)
reclassWrite <- function(type) {
  reclass = terra::classify(lr.crop, lands.meta[,c('category', type)])
  reproj = terra::project(reclass, 'EPSG:32735')
  terra::writeRaster(reproj,
                     filename=procpath(paste0(type, 'Raster.tif')),
                     overwrite=TRUE)
}
reclassWrite('movement')
reclassWrite('forage')
reclassWrite('shelter')

write.csv(lands.meta, file=rawpath('kaza_landcover', 'landcover_metadata.csv'))
