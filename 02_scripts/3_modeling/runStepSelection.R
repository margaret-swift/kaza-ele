# runStepSelection.R
# Created 13 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# amt vignette: 
#   https://cran.r-project.org/web/packages/amt/vignettes/p4_SSF.html
# interpreting SSF results: 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2656.13441


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/runStepSelection.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(sp, terra, amt)
setDataPaths('elephant')
load(procpath('eleKhau.rdata'))
setDataPaths('geographic')
load(procpath('geographic.rdata'))
load(procpath('stepSelectionParams.rdata'))


# ******************************************************************************
#                              RUN A SIMPLE SSF
# ******************************************************************************

m0 <- ssf.df |> fit_clogit(case_ ~ cover + strata(step_id_))
m1 <- ssf.df |> fit_clogit(case_ ~ cover + evi + strata(step_id_))
m2 <- ssf.df |> fit_clogit(case_ ~ cover + evi + waterdist + strata(step_id_))
summary(m0)
summary(m1)
summary(m2)





















# ******************************************************************************
#                               CREATE DATA 
# ******************************************************************************

# # Only run this section the first time
# khau.new <- st_transform(khau, "EPSG:32734") %>% as("Spatial")
# ssdf.khau <- ele.khau %>% dplyr::select(x, y, DATE.TIME) %>%
#   rename(x_=x, y_=y, t_=DATE.TIME) %>%
#   as_track()
# ssdf.khau$burst_ <- ele.khau$BURST

# # read in and crop rasters
# er <- rast(rawpath('khaudum_landsat_rasters', 'waterdist_evi_raster_dry_2017.tif'))
# lr <- rast(rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif'))
# names(lr) <- 'category'
# lr.khau <- crop(lr, khau.new)
# 
# # summarize sampling rate... is this right??
# summarize_sampling_rate(ssdf.khau)
# 
# # set up data for SSF by creating random steps and adding covariate values
# lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))
# ssf.df <- ssdf.khau %>%  
#   steps_by_burst() %>% 
#   random_steps(n_control = 15) %>% 
#   extract_covariates(lr) %>% 
#   extract_covariates(er) %>% 
#   left_join(lands.meta, by='category') %>% 
#   relocate(category, .before=class) %>% 
#   mutate(veg = factor(veg), cover = factor(cover)) %>% 
#   dplyr::select(-X, -sqkm) %>% 
#   rename_with(tolower)
# save(khau.new, ssdf.khau, ssf.df, file=procpath('stepSelectionParams.rdata'))
