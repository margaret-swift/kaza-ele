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

# load geographic data
setDataPaths('geographic')
load(procpath('geographic.rdata'))

# landcover types metadata
lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))
khau.new <- st_transform(khau, "EPSG:32734") %>% as("Spatial")

# load file with updated random_steps.bursted_steps_xyt; not sure why 
#   this one works better but w/e. added a progress bar too.
source(here::here('02_scripts','00_support', 'random_steps.R'))

# ******************************************************************************
#                                 RASTER DATA 
# ******************************************************************************

# hillslope raster
hr <- rast(rawpath('hillslope_khau.tif'))
names(hr) <- 'hillslope'

# EVI raster
er <- rast(rawpath('khaudum_landsat_rasters', 'waterdist_evi_raster_dry_2017.tif'))

# landcover data from WWF
lr <- rast(rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif'))
names(lr) <- 'category'

# set cover level order
cover_levels <- c('cropland', 'bare/water', 'open', 'sparse', 'closed')
veg_levels <- c('nonveg', 'cropland', 'bushland', 'herbaceous/wet', 'forest/woodland')


# ******************************************************************************
#                               CREATE DATA HSF
# ******************************************************************************

# create tracks from xy data
dat <- ele.khau %>% 
  arrange(INX) %>% 
  dplyr::select(INX, ID, x, y, DATE.TIME, TOD, BURST,
                SZN_2, SZN_4, SZN_6, STATE) %>%
  rename(x_=x, y_=y, t_=DATE.TIME)

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
#                               CREATE DATA SSF
# ******************************************************************************

ssf.df <- base_data %>% 
  # extract raster information
  extract_covariates(lr) %>%
  extract_covariates(er) %>%
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

# TODO
# remove alternate points that might cross a barrier

# save that s!
save(ssf.df, file=procpath('stepSelectionParamsHSF.rdata'))

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

