# cleanGeographic.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanGeographic.R')
source(here::here('02_scripts', 'utilities.R'))
setDataPaths('geographic')
load(here(procpath, 'geographic.rdata'))
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
     file=here(procpath, 'geographic.rdata'))

# EOF
