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

here::i_am('02_scripts/1_cleaning/createStepSelection.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(sp, terra, amt)
setDataPaths('elephant')
load(procpath('eleKhau.rdata'))
setDataPaths('geographic')
load(procpath('geographic.rdata'))
lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))

# ******************************************************************************
#                               RASTER DATA 
# ******************************************************************************

khau.new <- st_transform(khau, "EPSG:32734") %>% as("Spatial")
# read in and crop rasters
hr <- rast(rawpath('hillslope_khau.tif'))
er <- rast(rawpath('khaudum_landsat_rasters', 'waterdist_evi_raster_dry_2017.tif'))
lr <- rast(rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif'))
names(lr) <- 'category'
names(hr) <- 'hillslope'
lr.khau <- crop(lr, khau.new)

# set cover levels
coverlvl <- c('bare/water', 'cropland', 'open', 'sparse', 'closed')


# ******************************************************************************
#                               CREATE DATA 
# ******************************************************************************

# create tracks from xy data
ssdf.khau <- ele.khau %>% 
  dplyr::select(INX, BURST, x, y, DATE.TIME, TOD,
                SZN_2, SZN_4, SZN_6) %>%
  rename(x_=x, y_=y, t_=DATE.TIME) %>%
  as_track(all_cols=TRUE) %>% 
  track_resample(rate = hours(1), tolerance = minutes(10)) 

# summarize sampling rate and remove low-sample bursts
rates <- summarize_sampling_rate_many(ssdf.khau, cols="BURST")
burst.rm <- rates$BURST[which(rates$n < 20)]

# run with 150 alternate points and save
ssf.df <- ssdf.khau %>% 
  filter(!(BURST %in% burst.rm)) %>%
  steps_by_burst(keep_cols="start") %>%
  random_steps(n_control = 150) %>%
  extract_covariates(lr, where="both") %>%
  # extract_covariates(hr, where="both") %>%
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
  mutate(cover_start=factor(cover_start, levels=coverlvl),
         cover_end=factor(cover_end, levels=coverlvl))

# # set up water distance bins
# binMe <- function(v, bindf) {
#   inx <- which((v >= bindf$binmin) + (v < bindf$binmax) == 2)[1]
#   return(bindf$name[inx])
# }
# waterbins <- data.frame(name = c('0-500m', '500-1500m', '1500-2500m',
#                                  '2500-5000m', '5000-7500m', '7500-10,000m',
#                                  '>10,000m'),
#                         binmin = c(0, 500, 1500, 2500, 5000, 7500, -99),
#                         binmax = c(500, 1500, 2500, 5000, 7500, 11000, 0))
# binWater <- function(v) binMe(v, waterbins)
# bins <- unlist( lapply( ssf.df$waterdist, binWater ) )
# ssf.df$waterbin <-  factor(bins, levels=waterbins$name)
# inx.na = ssf.df$waterdist == -99
# ssf.df$waterdist[inx.na] <- NA

# add weights
ssf.df$w <- ifelse(ssf.df$case_, 1, 5000)

save(khau.new, ssdf.khau, ssf.df, file=procpath('stepSelectionParamsISSF.rdata'))




# ******************************************************************************
#                           DETERMINE AVAILABLE POINTS
# ******************************************************************************

# set up data for SSF by creating random steps and adding covariate values
# see Frieberg et al 2021 for details.

# # set up df to save values
# range <- seq(1, 200, by=10)
# NR = length(range)
# coefs.df <- data.frame(matrix(0, nrow=NR, ncol=3))
# names(coefs.df) <- c('hillslope', 'evi', 'waterdist')
# 
# # run model and extract coefs
# tic()
# for (i in 1:NR) {
#   message('i: ', i)
#   data <- ssdf.khau %>%
#     steps_by_burst() %>%
#     random_steps(n_control = i) %>%
#     extract_covariates(lr) %>%
#     extract_covariates(hr) %>%
#     extract_covariates(er) %>%
#     left_join(lands.meta, by='category') %>%
#     relocate(category, .before=class) %>%
#     mutate(veg = factor(veg), cover = factor(cover)) %>%
#     dplyr::select(-X, -sqkm) %>%
#     rename_with(tolower)
#   mod <- fit_clogit(data, case_ ~ hillslope + evi + waterdist + strata(step_id_))
#   coefs.df[i,] <- mod$model$coefficients
# }
# coefs.df$npoints = range
# toc()
# 
# 
# # when fitted coefficients level off, npoints is good
# p1 <- ggplot(data=coefs.df, mapping=aes(x=npoints, y=hillslope)) +
#   geom_point(color='purple') +
#   geom_smooth(color='purple', stat="smooth") + plot.theme
# 
# p2 <- ggplot(data=coefs.df, mapping=aes(x=npoints, y=evi)) +
#   geom_point(color='green') +
#   geom_smooth(color='green', stat="smooth") + plot.theme
# 
# p3 <- ggplot(coefs.df, aes(x=npoints, y=waterdist)) +
#   geom_point(color='cyan') +
#   geom_smooth(color='cyan', stat='smooth') + plot.theme
# p1 / p2 / p3
