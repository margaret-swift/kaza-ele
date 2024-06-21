# createHMM_SSF.R
# Created 18 June 2024
# Margaret Swift <margaret.swift@cornell.edu>

# HMM SSF package info
# devtools::install_github("NJKlappstein/hmmSSF")
# hmmSSF vignette: 
# https://github.com/NJKlappstein/hmmSSF/blob/main/vignettes/hmmSSF_introduction.pdf

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

# packages and HERE
pacman::p_load(here, hmmSSF)
i_am('02_scripts/3_modeling/hmm_ssf/createHMM_SSF.R')

# utilities
source(here('02_scripts', 'utilities.R'))

# loading data
setDataPaths('geographic')
load(procpath('geographic.rdata'))
load(procpath('stepSelectionParamsHSF.rdata'))
load(here(outdir, 'hmm', 'hmm_ssf.rdata'))
aoi = st_read(rawpath('kaza_aoi', 'kaza_aoi.shp'))

# ******************************************************************************
#                                PREP DATA FOR HMM
# ******************************************************************************

# load data
library(momentuHMM)
setDataPaths('elephant')
load(procpath('ele.rdata'))

# choose only females and 1-hour fix rates
data <- ele.df %>% nog() %>%
  dplyr::filter(SEX == "F", FIXRATE == "1 hour") %>%
  mutate(SEASON=factor(SZN_4),
         ID=BURST,
         DTM = ifelse(START.COUNT, NA, DTM)) %>%
  dplyr::select(INX, ID, BURST, SEX, SEASON, DATE.TIME, X, Y, DIST, DTM)

# remove first and last few fixes for each burst
bursts <- unique(data$ID)
inx.rm <- c()
for (i in 1:length(bursts)) {
  burst <- bursts[i]
  d <- data[data$ID == burst,]
}


# Crawl Wrap: Fill in missing data steps (85 tracks, t = 3m30s)
crwOut <-crawlWrap(obsData  = data,
                   timeStep = "hour",
                   Time.name= "DATE.TIME",
                   coord=c('X', 'Y'),
                   ncore=3,
                   fillCols=TRUE,
                   theta=c(6.855, -0.007)
)
crwOut$crwPredict$date <- as.POSIXct(crwOut$crwPredict$DATE.TIME)

# predicted vs. real data plot -- looks like it's just shifted a little
# to regularize the time steps
ggplot(data=crwOut$crwPredict[1:1000,]) +
  geom_point(aes(x=mu.x, y=mu.y, color=locType), size=2, alpha=0.85)

# Prep Data: add step lengths and turning angles again
pdata <- prepData( crwOut,
                   type="UTM",
                   coordNames=c("X", "Y"),
                   covNames=c('SEASON')
)

hmm.data <- pdata %>%
  mutate(Animal.ID = floor(BURST/1000),
         Loc.Type = ifelse(is.na(INX), "simulated", "real")) %>%
  dplyr::select(ID, Animal.ID, BURST, Loc.Type,
                date, SEASON, x, y, step, angle, se.mu.x, se.mu.y
                )

# ******************************************************************************
#                           CREATE USED/AVAILABLE DATA
# ******************************************************************************
 
# create tracks from xy data
track <- hmm.data %>%
  rename(time=date) %>% 
  dplyr::select(ID, x, y, time)

data <- get_controls(obs = track, n_controls = 50, distr = "gamma")
xydat <- data[,c('x', 'y')]

# ******************************************************************************
#                               ADD LANDSCAPE DATA
# ******************************************************************************

setDataPaths('geographic')

# landcover data from WWF
lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))
lr <- rast(rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif'))

# extract LC to data; this should take about 30sec.
data$category <- terra::extract(lr, xydat)$layer

# set cover level order
cover_levels <- c('cropland', 'bare/water', 'open', 'sparse', 'closed')
veg_levels <- c('nonveg', 'cropland', 'bushland', 'herbaceous/wet', 'forest/woodland')
hmm.ssf.dat <- data %>%
  # join classes from lands.meta
  rename_with(tolower) %>%
  left_join(lands.meta[,c('category', 'class')], by='category') %>%
  left_join(lands.meta[,c('category', 'cover_class.x')], by='category') %>%
  left_join(lands.meta[,c('category', 'veg_class.x')], by='category') %>%
  rename("cover_class"="cover_class.x",
         "veg_class"="veg_class.x") %>%
  # reorder factor levels
  mutate(cover_class=factor(cover_class, levels=cover_levels),
         veg_class = factor(veg_class, levels=veg_levels),
         category = factor(category), 
         class=factor(class))


# ******************************************************************************
#                                   EXTRACT EVI
# ******************************************************************************
setDataPaths('geographic')
hmm.ssf.dat$evi <- hmm.ssf.dat$season <- NA

# set up dry and wet season boundaries
makeSznInt <- function(y1, y2, d1, d2) { 
  dd1 = as.Date(paste(y1, d1, sep="-"))
  dd2 = as.Date(paste(y2, d2, sep="-"))
  message('  pulling data from: ', dd1, ' - ', dd2)
  return(interval(dd1, dd2))
}

# loop over rasters
slices <- c('10-01', '01-15', '05-01', '07-15')
szn.df <- data.frame(start=slices, end=slices[c(2:4, 1)])
row.names(szn.df) <- c('wet_1', 'wet_2', 'dry_1', 'dry_2')

rastnames <- list.files(rawpath('evi_rasters', 'evi_kaza_aoi'), full.names=TRUE)
data <- hmm.ssf.dat
for (f in rastnames) {
  evi_rast <- rast(f)
  fname <- gsub('.*evi_', '', f)
  message('running raster ', fname)
  
  # pull year and season to get date ranges
  y2 <- as.numeric( gsub('_.*', '', fname) )
  szn=gsub('.*[0-9]{4}_|.tif$', '', fname)
  y1 = ifelse(szn=='wet_1', y2-1, y2)
  ds <- szn.df[szn,]
  my.inx <- data$time %within% makeSznInt(y1,y2,ds$start,ds$end)
  
  # extract covariates and save
  evi.dat <- terra::extract(evi_rast, xydat[my.inx,])$EVI
  message('  checksum: ', length(evi.dat) == sum(my.inx))
  data$evi[my.inx] <- evi.dat
  data$season[my.inx] <- szn
}
data <- data %>% 
  # while we're working with landsat 8 we have to be mindful of availability
  filter(time > as.Date('2013-05-01'),
         time < as.Date('2022-10-01')) 

# make sure all points are within the dataset
crs='EPSG:32734'
aoi_projected = st_transform(aoi, crs=crs)
dat_projected = data %>% 
  sample_n(50000) %>% 
  st_as_sf(coords=c('x', 'y'), crs='EPSG:32734') %>% 
  st_intersection(aoi_projected) %>% 
  mutate(YEAR = as.factor(year(time)))
ggplot() + 
  geom_sf(data=aoi_projected, fill=NA, color='black') +
  geom_sf(data=dat_projected, 
          mapping=aes(color=evi)) + 
  scale_color_distiller(palette="Greens", direction=1) + 
  facet_wrap(~YEAR)

hmm.ssf.dat <- data

# ******************************************************************************
#                               SAVE DATA FILES
# ******************************************************************************
save(hmm.ssf.dat, file=here(outdir, "hmm", "hmm_ssf.rdata"))



# EOF