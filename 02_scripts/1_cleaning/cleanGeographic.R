# cleanGeographic.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanGeographic.R')
source(here::here('02_scripts', 'utilities.R'))
setDataPaths('geographic')
load(procpath('geographic.rdata'))
my_crs <- 4326

loadData <- function(path) {
  file <- list.files(here(rawpath, path), pattern='.shp$', full.names=TRUE)
  data <- st_read(file) %>% st_transform(crs=my_crs)
  data
}

# Pulling all data
fences <- loadData('merged_fences')
roads <- loadData('primary_trunk_rds_5countries')
rivers <- loadData('rivers_main_kaza_dig')
kaza <- loadData('kaza_boundary')
khau <- loadData('khaudum_boundary')
waters_art <- loadData('khaudum_artificial_waters')
waters_nat <- loadData('khaudum_natural_waters')

save(fences, roads, rivers,
     kaza, khau, waters_art, waters_nat, 
     file=procpath('geographic.rdata'))


# ******************************************************************************
#                        Set up raster resistance layer
# ******************************************************************************

load(procpath('landcover_rasters.rdata'))

# kinda the only way I could find to mosaic rasters with different extents...
# https://stackoverflow.com/questions/67169266/error-in-mosaic-of-rasters-from-different-extent-using-terra-package-in-r
f1 <- here('03_output', 'tmp', 'test1.tif')
f2 <- here('03_output', 'tmp', 'test2.tif')
r1 <- writeRaster(lands.raster.1, f1, overwrite=TRUE)
r2 <- writeRaster(lands.raster.2, f2, overwrite=TRUE)
lands.vrt <- vrt(c(f1, f2), "03_output/tmp/test.vrt", overwrite=TRUE) 

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

forage.mat <- classify(lands.vrt, lands.meta[,c('category', 'forage')]) %>%
  as.matrix(wide=TRUE)
move.mat <- classify(lands.vrt, lands.meta[,c('category', 'movement')]) %>%
  as.matrix(wide=TRUE)
shelter.mat <- classify(lands.vrt, lands.meta[,c('category', 'shelter')]) %>%
  as.matrix(wide=TRUE)

save(shelter.mat, move.mat, forage.mat, file=procpath('landsmats.rdata'))

