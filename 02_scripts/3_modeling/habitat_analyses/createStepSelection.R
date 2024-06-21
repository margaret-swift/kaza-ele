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

# load ele data
setDataPaths('elephant')
load(procpath('eleKhau.rdata'))
load(here(outdir, 'hmm', 'hmm_ssf.rdata'))


# load geographic data
setDataPaths('geographic')
load(procpath('geographic.rdata'))

# landcover types metadata
lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))
khau.new <- st_transform(khau, "EPSG:32734") %>% as("Spatial")
khau.tall.new <-  st_transform(khau_tall, "EPSG:32734")

# load file with updated random_steps.bursted_steps_xyt; not sure why 
#   this one works better but w/e. added a progress bar too.
source(here::here('02_scripts','00_support', 'random_steps.R'))

# ******************************************************************************
#                                 RASTER DATA 
# ******************************************************************************

# hillslope raster
# hr <- rast(rawpath('hillslope_khau.tif'))
# names(hr) <- 'hillslope'

# water distance raster
wdr <- rast(rawpath('khaudum_landsat_rasters', 
                   'waterdist_evi_raster_dry_2017.tif'))

# landcover data from WWF
lr <- rast(rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif'))
names(lr) <- 'category'

# set cover level order
cover_levels <- c('cropland', 'bare/water', 'open', 'sparse', 'closed')
veg_levels <- c('nonveg', 'cropland', 'bushland', 'herbaceous/wet', 'forest/woodland')

# filter ele data to just that of the tall AOI
x <- ele.khau %>% 
  filter(!is.na(x), !is.na(y)) %>% 
  st_as_sf(coords=c('x', 'y'), crs="EPSG:32734")
ele.khau.tall <- st_intersection(x, khau.tall.new)

# ggplot() + 
#   geom_sf(data=khau.tall.new) +
#   geom_sf(data=ele.khau.tall)

# ******************************************************************************
#                           CREATE USED/AVAILABLE DATA
# ******************************************************************************

# create tracks from xy data
dat <- ele.khau.tall %>% 
  arrange(INX) %>% 
  dplyr::select(INX, ID, X, Y, DATE.TIME, TOD, BURST,
                SZN_2, SZN_4, SZN_6, STATE) %>%
  rename(x_=X, y_=Y, t_=DATE.TIME)

bursts <- unique(dat$BURST)
ssdf.khau <- NULL
for (i in 1:length(bursts)) {
  b <- bursts[i]
  xyt <- dat %>% filter(BURST == b)
  track <- as_track(xyt, all_cols=TRUE) %>% 
    track_resample(rate = hours(1), tolerance = minutes(10)) 
  if (is.null(ssdf.khau)) ssdf.khau <- track
  else ssdf.khau <- rbind(ssdf.khau, track)
}

# summarize sampling rate and remove low-sample bursts
rates <- ssdf.khau %>% group_by(burst_) %>% summarize(n=n())
burst.rm <- rates$burst_[which(rates$n < 20)]
ssdf.khau <-  filter(ssdf.khau, !(burst_ %in% burst.rm))

# For some reason, INX 505934 is a problem child (sl is SO large). 
#  But, if I remove it before running steps-by-burst, 505933 becomes problem child.
pc = 505934
n_control = 150
w = 5000

base.data <- ssdf.khau %>%
  steps_by_burst(keep_cols="start") %>%
  filter(INX != pc) %>% 
  # run with 150 alternate points
  random_steps.bursted_steps_xyt(n_control = n_control) %>%
  # reweight alternate points
  mutate(w = ifelse(case_, 1, w))



# ******************************************************************************
#                               EXTRACT LANDCOVER
# ******************************************************************************

ssf <- base.data %>% 
  # extract raster information
  extract_covariates(lr) %>%
  # join classes from lands.meta
  rename_with(tolower) %>% 
  left_join(lands.meta[,c('category', 'cover_class.x')], by='category') %>% 
  left_join(lands.meta[,c('category', 'veg_class.x')], by='category') %>% 
  rename("cover_class"="cover_class.x",
         "veg_class"="veg_class.x") %>% 
  # reorder factor levels
  mutate(cover_class=factor(cover_class, levels=cover_levels),
         veg_class = factor(veg_class, levels=veg_levels),
         state = factor(ifelse(state == "correlated walk", "exploring", state),
                         levels=c('resting', 'foraging', 'exploring')),
         category = factor(category))

ssf$season <- ssf$settle_dist <- ssf$evi <- NA

# ******************************************************************************
#                           EXTRACT WATER DISTANCE
# ******************************************************************************

# pull raster and extract covariates
ssf <- extract_covariates(ssf, wdr$waterDist) %>% rename(water_dist = waterDist)

# ******************************************************************************
#                                 EXTRACT EVI
# ******************************************************************************
# set up dry and wet season boundaries
dry.int <- c('10-01', '05-01')
wet.int <- c('05-01', '09-30')

# loop over rasters
rastnames <- list.files(rawpath('evi_rasters'), pattern=".tif$", full.names=TRUE)
for (f in rastnames) {
  evi_rast <- rast(f)
  fname <- gsub('.*evi_', '', f)
  message('running raster ', fname)
  
  # pull year and season to get date ranges
  yr <- as.numeric( gsub('_.*', '', fname) )
  szn <- gsub('.*_|.tif', '', fname)
  if (szn == "dry") {
    d1 = paste(yr-1, dry.int[1], sep="-")
    d2 = paste(yr, dry.int[2], sep="-")
  } else {
    d1 = paste(yr, wet.int[1], sep="-")
    d2 = paste(yr, wet.int[2], sep="-")
  }
  my.inx <- ssf$t1_ %within% interval(as.Date(d1), as.Date(d2))
  
  # extract covariates and save
  evi.vals <- extract_covariates(ssf[my.inx,], evi_rast)
  ssf$evi[my.inx] <- evi.vals$EVI
  ssf$season[my.inx] <- szn
}
ssf <- ssf %>% filter(!is.na(evi)) 


# ******************************************************************************
#                             EXTRACT BURNED AREA
# ******************************************************************************


## TODO: 
## -- figure out how to handle burn dates ('burned in last month'?)
## -- transform MODIS data to UTM smh :( 

rastnames <- list.files(rawpath('MODIS Burned area'), pattern=".tif$", full.names=TRUE)
ssf$DOY = yday(ssf$t1_)
ssf$was_burned <- FALSE

for (i in 1:length(rastnames)) {
  f <- rastnames[i]
  modis_rast <- rast(f) %>% terra::project("EPSG:32734")
  fname <- gsub('.*export_|.tif$', '', f)
  message('running raster ', fname)
  
  # pull year and season to get date ranges
  yr <- as.numeric( gsub('.*export_fire_|.tif', '', f) )
  d1 = as.Date(paste(yr, "01-01", sep="-"))
  d2 = as.Date(paste(yr+1, "01-31", sep="-")) # allow for 30-day window into next year
  my.inx <- ssf$t1_ %within% interval(d1, d2)
  
  # extract covariates and save
  if (any(my.inx)) {
    modis.vals <- extract_covariates(ssf[my.inx,], modis_rast) %>% 
      mutate(DOY = ifelse(year(t1_) == yr+1, DOY + 365, DOY))
    inx.rm <- modis.vals$BurnDate == 0
    was_burned <- round(modis.vals$DOY - modis.vals$BurnDate) %in% 1:30
    was_burned[inx.rm] <- FALSE
    ssf$was_burned[my.inx] <- ssf$was_burned[my.inx] + was_burned
  }
}
ssf$was_burned <- ssf$was_burned > 0
ssf <- ssf %>% filter(year(t1_) < 2021)



# ******************************************************************************
#                               SAVE SSF DATA
# ******************************************************************************


# TODO
# remove alternate points that might cross a barrier

# save that s!
save(ssf, file=procpath('stepSelectionParamsHSF.rdata'))

# ******************************************************************************
#                               CREATE DATA ISSF
# ******************************************************************************

issf.df <- base_data %>% 
  extract_covariates(lr, where="both") %>%
  extract_covariates(er, where="both") %>%
  rename_with(tolower) %>% 
  # match to start category
  rename(category=category_start) %>% 
  left_join(lands.meta[,c('category', 'cover_class')], by='category') %>% 
  rename(cover_start=cover_class,
         category_start=category,
         category=category_end) %>% 
  left_join(lands.meta[,c('category', 'cover_class')], by='category') %>% 
  rename(cover_end=cover_class,
         category_end=category) %>% 
  relocate(cover_start, .after=category_start) %>%
  relocate(cover_end, .after=category_end) %>%
  mutate(cover_start=factor(cover_start, levels=cover_levels),
         cover_end=factor(cover_end, levels=cover_levels))

save(issf.df, file=procpath('stepSelectionParamsISSF.rdata'))

