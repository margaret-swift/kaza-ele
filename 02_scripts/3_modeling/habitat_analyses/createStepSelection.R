# createStepSelection.R
# Created 13 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# amt vignette: 
#   https://cran.r-project.org/web/packages/amt/vignettes/p4_SSF.html
# interpreting SSF results: 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2656.13441


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/habitat_analyses/createStepSelection.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(sp, terra, amt)

# load linear, boundary, and elephant data
quickload()
quickload('elephant')

# load file with updated random_steps.bursted_steps_xyt; not sure why 
#   this one works better but w/e. added a progress bar too.
source(here::here('02_scripts','00_support', 'random_steps.R'))

# loading existing file if need be
# setOutPath('habitat_selection')
# load(outpath('stepSelectionParamsHSF_DF.rdata'))


# ******************************************************************************
#                                 RASTER DATA 
# ******************************************************************************

# water distance raster
setDataPaths('waters')
wdr <- rast(rawpath('waterdist_evi_raster_dry_2017',
                   'waterdist_evi_raster_dry_2017.tif'))

# landcover data from WWF
setDataPaths('landcover')
lands.meta <- read.csv(metapath('landcover_meta_2005.csv'))
lr <- raster::raster(rawpath('landcover_2005_fixed.tif'))

# set cover level order
cover_levels <- c('nonveg', 'open', 'closed', 'sparse')
veg_levels <- c('nonveg', 'wetland', 'cropland', 'grassland', 'woodland', 'bushland')
struc_levels <- c('nonveg', 'grassland','cropland',
                  'wetland', 'omiramba', 
                  'open bush', 'sparse bush', 'closed bush')

# ******************************************************************************
#                           CREATE USED/AVAILABLE DATA
# ******************************************************************************


# create tracks from xy data
dat <- ele.df %>% 
  arrange(INX) %>% 
  filter(SEX=="F") %>% 
  dplyr::select(INX, ID, X, Y, DATE.TIME, TOD, BURST,
                SZN_2, SZN_4, SZN_6) %>%
  rename(x_=X, y_=Y, t_=DATE.TIME)
# x = sample_n(dat, 10000)
# plot(x$x_, x$y_)

bursts <- unique(dat$BURST)
ssdf <- NULL
for (i in 1:length(bursts)) {
  b <- bursts[i]
  message(b)
  xyt <- dat %>% filter(BURST == b)
  track <- as_track(xyt, all_cols=TRUE) %>% 
    track_resample(rate = hours(1), tolerance = minutes(10)) 
  track$burst_ = (i*1000) + track$burst_
  if (is.null(ssdf)) ssdf <- track
  else ssdf <- rbind(ssdf, track)
}

# summarize sampling rate and remove low-sample bursts
rates <- ssdf %>% group_by(burst_) %>% summarize(n=n())
burst.rm <- rates$burst_[which(rates$n < 20)]
ssdf <-  filter(ssdf, !(burst_ %in% burst.rm))

# set number of control steps and weights
nc = 100
w = 5000
base.data <- ssdf %>%
  steps_by_burst(keep_cols="start") %>%
  # add n_control alternate points
  .random_steps.bursted_steps_xyt(n_control=nc) %>%
  # reweight alternate points
  mutate(w = ifelse(case_, 1, w))

steps.na <- base.data %>% 
  filter(is.na(x2_)) %>% 
  dplyr::select(step_id_) %>% 
  unique()
base.data$FLAGME <- base.data$step_id_ %in% steps.na$step_id_


# ******************************************************************************
#                               EXTRACT LANDCOVER
# ******************************************************************************

# Get xy data locations in Alberts 1984 projection to match WWF landcover data
xy.dat <- base.data %>% nog() %>% 
  filter(!FLAGME) %>% 
  dplyr::select('x2_', 'y2_') %>% 
  as.data.frame()
names(xy.dat) <- c("X", "Y")
xy.dat.utm <- st_as_sf(xy.dat, coords=c('X', 'Y'))
st_crs(xy.dat.utm) <- "EPSG:32734"
xy.dat.alb <- st_transform(xy.dat.utm, crs="ESRI:102022")

# extract raster information
ssf.lc <- xy.dat.alb
x <- terra::extract(lr, xy.dat.alb)
ssf.lc$category = x
  # join classes from lands.meta
ssf.lc <- ssf.lc %>% 
  nog() %>% 
  left_join(lands.meta, by='category') %>% 
  mutate(cover=factor(cover, levels=cover_levels),
         vegetation = factor(vegetation, levels=veg_levels),
         structure = factor(structure, levels=struc_levels),
         category = factor(category)) %>% 
    dplyr::select(cover, structure, vegetation)

# ******************************************************************************
#                           EXTRACT WATER DISTANCE
# ******************************************************************************

# pull raster and extract covariates
setDataPaths('waters')
ssf.wat <- terra::extract(wdr$waterDist, xy.dat.utm) %>% 
  rename(water_dist = waterDist)


# ******************************************************************************
#                                 EXTRACT EVI
# ******************************************************************************

# set up dry and wet season boundaries
setDataPaths('evi')
slices <- c('10-01', '01-15', '05-01', '07-15')
szn.df <- data.frame(start=slices, end=slices[c(2:4, 1)])
row.names(szn.df) <- c('wet_1', 'wet_2', 'dry_1', 'dry_2')
makeSznInt <- function(y1, y2, d1, d2) { 
  interval( as.Date(paste(y1, d1, sep="-")), 
            as.Date(paste(y2, d2, sep="-"))) 
}

# loop over rasters
rastnames <- list.files(rawpath('evi_kaza_aoi', 'evi_kaza_aoi'), 
                        full.names=TRUE)
ssf.evi <- base.data %>% filter(!FLAGME)
for (f in rastnames) {
  evi_rast <- rast(f)
  fname <- gsub('.*evi_', '', f)
  message('running raster ', fname)
  
  # pull year and season to get date ranges
  y2 <- as.numeric( gsub('_.*', '', fname) )
  szn=gsub('.*[0-9]{4}_|.tif$', '', fname)
  y1 = ifelse(szn=='wet_1', y2-1, y2)
  ds <- szn.df[szn,]
  my.inx <- ssf.evi$t1_ %within% makeSznInt(y1,y2,ds$start,ds$end)
  
  # extract covariates and save
  evi.dat <- terra::extract(evi_rast, xy.dat.utm[my.inx,])$EVI
  message('  checksum: ', length(evi.dat) == sum(my.inx))
  ssf.evi$evi[my.inx] <- evi.dat
  ssf.evi$season[my.inx] <- szn
}


# ******************************************************************************
#                             EXTRACT BURNED AREA
# ******************************************************************************


## TODO: 
## -- figure out how to handle burn dates ('burned in last month'?)
## -- transform MODIS data to UTM smh :( 

# rastnames <- list.files(rawpath('MODIS Burned area'), pattern=".tif$", full.names=TRUE)
# ssf$DOY = yday(ssf$t1_)
# 
# for (i in 1:length(rastnames)) {
#   f <- rastnames[i]
#   modis_rast <- rast(f) %>% terra::project("EPSG:32734")
#   fname <- gsub('.*export_|.tif$', '', f)
#   message('running raster ', fname)
#   
#   # pull year and season to get date ranges
#   yr <- as.numeric( gsub('.*export_fire_|.tif', '', f) )
#   d1 = as.Date(paste(yr, "01-01", sep="-"))
#   d2 = as.Date(paste(yr+1, "01-31", sep="-")) # allow for 30-day window into next year
#   my.inx <- ssf$t1_ %within% interval(d1, d2)
#   
#   # extract covariates and save
#   if (any(my.inx)) {
#     modis.vals <- extract_covariates(ssf[my.inx,], modis_rast) %>% 
#       mutate(DOY = ifelse(year(t1_) == yr+1, DOY + 365, DOY))
#     inx.rm <- modis.vals$BurnDate == 0
#     was_burned <- round(modis.vals$DOY - modis.vals$BurnDate) %in% 1:30
#     was_burned[inx.rm] <- FALSE
#     ssf$was_burned[my.inx] <- ssf$was_burned[my.inx] + was_burned
#   }
# }
# ssf$was_burned <- ssf$was_burned > 0
# ssf <- ssf %>% filter(year(t1_) < 2021)
# 

# ******************************************************************************
#                               SAVE SSF DATA
# ******************************************************************************

# TODO
# remove alternate points that might cross a barrier

ssf <- cbind(ssf.evi, ssf.lc, ssf.wat) %>% 
  nog() %>% dplyr::select(-geometry) %>% 
  # while we're working with landsat 8 we have to be mindful of availability
  filter(t1_ > as.Date('2013-05-01'),
         t1_ < as.Date('2022-10-01')) %>% 
  mutate(cos_ta = cos(ta_),
         log_sl = log(sl_))

# save that s!
setOutPath('habitat_selection')
save(ssf, file=outpath('stepSelectionParamsHSF_DF.rdata'))
