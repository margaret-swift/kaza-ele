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
pacman::p_load(here, hmmSSF, momentuHMM, 
               daiR, googleCloudStorageR)
i_am('02_scripts/3_modeling/hmm_ssf/createHMM_SSF.R')

# utilities
source(here('02_scripts', 'utilities.R'))

# loading data
# quickload()
quickload('elephant')
# setOutPath('habitat_selection')
# load(outpath('ssf_data', 'stepSelectionParamsHSF_DF.rdata'))

# loading data if you've already created it
setOutPath('hmm')
outpath('hmm_ssf.rdata')

# load(here(outdir, 'hmm', 'hmm_ssf.rdata'))
setDataPaths('boundaries')
aoi = st_read(rawpath('aoi_bots_nam', 'aoi.shp'))

# ******************************************************************************
#                                PREP DATA FOR HMM
# ******************************************************************************

# choose only females and 1-hour fix rates
data <- ele.df %>% nog() %>%
  dplyr::filter(SEX == "F", FIXRATE == "1 hour") %>%
  mutate(SEASON=factor(SZN_4),
         ID=BURST,
         DTM = ifelse(START.COUNT, NA, DTM)) %>%
  dplyr::select(INX, ID, BURST, SEX, SEASON, DATE.TIME, X, Y, DIST, DTM)

# remove first and last few fixes for each burst
bursts.keep <- data.frame(table(data$ID, dnn=c('BURST'))) %>% 
  filter(Freq > 4) %>% 
  dplyr::select(BURST)

# Crawl Wrap: Fill in missing data steps (85 tracks, t = 3m30s)
# Crawl Wrap: Fill in missing data steps (100 tracks, t = 4m40s)
freemem()
tic()
crwOut <-crawlWrap(obsData  = data %>% filter(ID %in% bursts.keep$BURST),
                   timeStep = "hour",
                   Time.name= "DATE.TIME",
                   coord=c('X', 'Y'),
                   ncore=3,
                   fillCols=TRUE,
                   theta=c(6.855, -0.007) #this is from the tutorial, not sure how to choose
)
toc()
rm(data)
freemem(TRUE)
gc()
crwOut$crwPredict$date <- as.POSIXct(crwOut$crwPredict$DATE.TIME)

# predicted vs. real data plot -- looks like it's just shifted a little
# to regularize the time steps
ggplot(data=crwOut$crwPredict[1:1000,]) +
  geom_point(aes(x=mu.x, y=mu.y, color=locType), size=2, alpha=0.85)
setOutPath('hmm')
save(crwOut, file=outpath('crwout.rdata'))

# Prep Data: add step lengths and turning angles again
tic()
pdata <- prepData( crwOut,
                   type="UTM",
                   coordNames=c("X", "Y"),
                   covNames=c('SEASON')
)

hmm.data <- pdata %>%
  mutate(Animal.ID = floor(BURST/1000),
         Loc.Type = ifelse(is.na(INX), "simulated", "real")) %>%
  dplyr::select(ID, Animal.ID, BURST, Loc.Type,
                date, SEASON, x, y, step, angle, se.mu.x, se.mu.y)
toc()

rm(crwOut)
rm(pdata)
setOutPath('hmm')
save(hmm.data, file=outpath('hmm_ssf.rdata'))




# ******************************************************************************
#                           CREATE USED/AVAILABLE DATA
# ******************************************************************************

tic()
track <- hmm.data %>%
  rename(time=date) %>% 
  dplyr::select(ID, x, y, time)
ssf.track <- get_controls(obs = track, n_controls = 50, distr = "gamma")
toc()
save(ssf.track, file=outpath('hmm_ssf_track.rdata'))


# ******************************************************************************
#                             CREATE XY PROJECTIONS
# ******************************************************************************

tic()
xy.dat <- ssf.track[,c('x', 'y')]
xy.dat.utm <- st_as_sf(xy.dat, coords=c('x', 'y'))
st_crs(xy.dat.utm) <- "EPSG:32734"
xy.dat.alb <- st_transform(xy.dat.utm, crs="ESRI:102022")
toc()
save(xy.dat.alb, file=outpath('hmm_ssf_xydat.rdata'))

## LOAD THESE
# load(outpath('hmm_ssf_track.rdata'))
# load(outpath('hmm_ssf_xydat.rdata'))

# ******************************************************************************
#                               ADD LANDSCAPE DATA
# ******************************************************************************

# set up proper projection data
ssf.lc <- xy.dat.alb

# landcover data from WWF (2005, not 2021)
setDataPaths('landcover')
lands.meta <- read.csv(metapath('landcover_meta_2005.csv'))
lr <- raster::raster(rawpath('kaza_landcover_2005', 'landcover_2005_fixed.tif'))

# set cover level order
cover_levels <- c('nonveg', 'open', 'closed', 'sparse')
veg_levels <- c('nonveg', 'wetland', 'cropland', 'grassland', 'woodland', 'bushland')
struc_levels <- c('nonveg', 'grassland','cropland',
                  'wetland', 'omiramba', 
                  'open bush', 'sparse bush', 'closed bush')

# extract raster information
lr.dat <- terra::extract(lr, ssf.lc)
ssf.lc$category = lr.dat

# join classes from lands.meta
ssf.lc <- ssf.lc %>% 
  nog() %>% 
  left_join(lands.meta, by='category') %>% 
  mutate(cover =      factor(cover,     levels=cover_levels),
         vegetation = factor(vegetation,levels=veg_levels),
         structure =  factor(structure, levels=struc_levels),
         category =   factor(category),
         class =      factor(class))


# ******************************************************************************
#                                 EXTRACT EVI
# ******************************************************************************

# set up proper projection data
ssf.evi <- cbind(ssf.track$time, xy.dat.utm)
ssf.evi$evi <- ssf.evi$season <- NA
names(ssf.evi)[1] <- 'time'

# load data
src <- "local" #"gcs" for cloud storage, otherwise "local"
if (src == "local") {
  setDataPaths('evi')
  rastnames <- list.files(rawpath('evi_okavango'), full.names=TRUE)
  
} else if (src == "gcs") {
  # load data from GCS
  rastnames <- gcs_list_objects(bucket="kaza_evi_rasters")[,1]
}


# set up dry and wet season boundaries
sl <- c('11-01', '01-01', '03-01', '05-01', '07-01', '09-01')
szn.df <- data.frame(start=sl, end=sl[c(2:length(sl), 1)])
row.names(szn.df) <- c('wet_1', 'wet_2', 'wet_3', 'dry_1', 'dry_2', 'dry_3')
makeSznInt <- function(y1, y2, d1, d2) { 
  interval( as.Date(paste(y1, d1, sep="-")), 
            as.Date(paste(y2, d2, sep="-"))) 
}

# loop over rasters
for (f in rastnames) {
  if (src == "local") {
    evi_rast <- terra::rast(f)
  } else if (src == "gcs") {
    # download file as tmp and open it that way
    # 
  }
  
  
  fname <- gsub('.*evi_', '', f)
  message('running raster ', fname)
  
  # pull year and season to get date ranges
  y2 <- as.numeric( gsub('_.*', '', fname) )
  szn=gsub('.*[0-9]{4}_|.tif$', '', fname)
  y1 = ifelse(szn=='wet_1', y2-1, y2)
  ds <- szn.df[szn,]
  my.inx <- ssf.evi$time %within% makeSznInt(y1,y2,ds$start,ds$end)
  
  # extract covariates and save
  evi.dat <- terra::extract(evi_rast, xy.dat.utm[my.inx,])$EVI
  message('  checksum: ', length(evi.dat) == sum(my.inx))
  ssf.evi$evi[my.inx] <- evi.dat
  ssf.evi$season[my.inx] <- szn
}

inx.rm <- (ssf.evi$evi < 0) + (ssf.evi$evi > 1.0) > 0
ssf.evi$evi[inx.rm] <- NA
ssf.evi <- ssf.evi %>% dplyr::select(-time, -season)
hist(ssf.evi$evi[!inx.rm], breaks=100)


# ******************************************************************************
#                              LAST CHECK
# ******************************************************************************

# while we're working with landsat 8 we have to be mindful of availability
hmm.ssf <- cbind(ssf.track, ssf.lc, ssf.evi) %>% 
  nog() %>% 
  filter(time > as.Date('2014-01-01'),
         time < as.Date('2022-10-01')) %>% 
  mutate(cos_ta = cos(angle),
         log_sl = log(step)) %>% 
  rename(id=ID)


# # make sure all points are within the dataset
# crs='EPSG:32734'
# aoi_projected = st_transform(aoi, crs=crs)
# dat_projected = hmm.ssf.dat %>% 
#   sample_n(50000) %>% 
#   st_as_sf(coords=c('x', 'y'), crs='EPSG:32734') %>% 
#   st_intersection(aoi_projected) %>% 
#   mutate(YEAR = as.factor(year(time)))
# ggplot() + 
#   geom_sf(data=aoi_projected, fill=NA, color='black') +
#   geom_sf(data=dat_projected, 
#           mapping=aes(color=evi)) + 
#   scale_color_distiller(palette="Greens", direction=1) + 
#   facet_wrap(~YEAR)

# ******************************************************************************
#                               SAVE DATA FILES
# ******************************************************************************
setOutPath('hmm')
save(hmm.ssf, file=outpath("hmm_ssf.rdata"))


# EOF