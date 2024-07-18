# addCovariates.R
# Created 04 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/addCovariates.R')
source(here::here('02_scripts', 'utilities.R'))
pacman::p_load(suncalc, zoo, lwgeom)

quickload()
quickload('elephant')
quickload('precipitation')
setDataPaths('precipitation')
# load(procpath('season_slices.rdata'))

# landcover data from WWF
setDataPaths('landcover')
lands.meta <- read.csv(metapath('landcover_meta_2005.csv'))
lr <- raster::raster(rawpath('landcover_2005_fixed.tif'))

# ******************************************************************************
#                               Precipitation
# ******************************************************************************

# add precip in the last week
precip.df$MM.RAIN_7 <- rollsum(precip.df$MM_RAIN, k=7, fill=NA, align="right")

p.inx <- match(ele.df$DATE, precip.df$DATE)
ele.df$MM.RAIN <- precip.df$MM_RAIN[p.inx]
ele.df$MM.RAIN_C <- precip.df$MM_RAIN_C[p.inx]
ele.df$MM.RAIN_7 <- precip.df$MM.RAIN_7[p.inx]


# ******************************************************************************
#                               Landcover type (WWF 2021)
# ******************************************************************************
# 
# # find landcover classes
# findLC <- function(val) {
#   inx <- ifelse(is.na(val), NA, 
#                 which.min(abs(lands.meta$category - val)))
#   return( inx )
# }
# 
# # this takes 19 seconds
# LC.df <- data.frame(INX=ele.df$INX, 
#                     lc1=raster::extract(lands.raster.1, ele.df),
#                     lc2=raster::extract(lands.raster.2, ele.df)) %>% 
#   mutate(LC_CLASS = ifelse(is.na(lc1), lc2, lc1))
# LC.df$LC_CATEG<- lands.meta$class[ unlist( lapply(LC.df$LC, function(e) findLC(e)) ) ]
# ele.df$LC_CLASS <- LC.df$LC_CLASS
# ele.df$LC_CATEG <- factor( LC.df$LC_CATEG, levels=lands.meta$class )

# ******************************************************************************
#                               Landcover type (WWF 2005)
# ******************************************************************************

# Get xy data locations in Alberts 1984 projection to match WWF landcover data
xy.dat <- ele.df %>% nog() %>% 
  dplyr::select('X', 'Y') %>% 
  as.data.frame()
xy.dat.utm <- st_as_sf(xy.dat, coords=c('X', 'Y'))
st_crs(xy.dat.utm) <- "EPSG:32734"
xy.dat.alb <- st_transform(xy.dat.utm, crs="ESRI:102022")

# extract raster information
lc <- xy.dat.alb
# x <- terra::extract(lr, xy.dat.alb)
lc$category = x
# join classes from lands.meta
lc <- lc %>% 
  nog() %>% 
  left_join(lands.meta, by='category') %>% 
  mutate(COVER=factor(cover, levels=cover_levels),
         VEGETATION = factor(vegetation, levels=veg_levels),
         STRUCTURE = factor(structure, levels=struc_levels),
         CATEGORY = factor(category),
         CLASS = factor(class)) %>% 
  dplyr::select(CATEGORY, CLASS, COVER, STRUCTURE, VEGETATION)
ele.df <- cbind(ele.df, lc)

# ******************************************************************************
#                             Adding time of day
# ******************************************************************************

data <- ele.df %>% nog() %>%
  rename(lat=LAT, lon=LON, date.time=DATE.TIME) %>%
  mutate(date = as.Date(as.POSIXct(date.time)))
sun.times <- getSunlightTimes(data=data[,c('date', 'lat', 'lon')])
findDiff <- function(e) difftime( e, data$date.time )

dawn <- difftime (sun.times$dawn, data$date.time ) %>% as.n()
sunriseEnd <- difftime (sun.times$sunriseEnd, data$date.time ) %>% as.n()
sunset <- difftime (sun.times$sunset, data$date.time ) %>% as.n()
night <- difftime (sun.times$night, data$date.time ) %>% as.n()
df <- data.frame(dawn, sunriseEnd, sunset, night)

# account for before dawn
signs <- sign(df) == -1
sums <- rowSums(signs)
sums[sums == 0] <- 4
time.names <- c('DAWN', 'DAY', 'DUSK', "NIGHT")
TOD <- time.names[sums]
ele.df$TOD <- TOD

# ******************************************************************************
#                            ASSIGNING SEASONS TO DATES
# ******************************************************************************

szn_date <- ele.df$DATE %>% 
  gsub('\\d{4}', '1899', .) %>% 
  gsub('02-29', '02-28', .) %>% # take care of leap years
  as.Date(.)
szn_week <- week(szn_date)
sznInt <- function(d1, d2, y='1899') {
  interval(as.Date(paste(y, d1, sep="-")), 
           as.Date(paste(y, d2, sep="-")))
}
# seasonal assignment - dry, wet
ele.df$SZN_2 <- ifelse(szn_date %within% sznInt("05-01", "10-31"), 'DRY', 'WET')

# seasonal shifts - wet, 
levels=c('WET1', 'WET2', 'DRY')
ele.df$SZN_3 <- factor(
  case_when(
    szn_week %in% 1:17 ~ levels[1],
    szn_week %in% 43:53 ~ levels[2],
    .default = levels[3]
  ), levels=levels)

# seasonal shifts - wet, dry, and two transitions
levels=c('T1', 'WET', 'T2', 'DRY')
ele.df$SZN_4 <- factor(
  case_when(
    szn_week %in% 43:48 ~ levels[1],
    szn_week %in% 12:17 ~ levels[3],
    szn_week %in% 18:42 ~ levels[4],
    .default = levels[2]
  ), levels=levels)

# ******************************************************************************
#                             Adding distance from water
# ******************************************************************************
# x = st_distance(ele.df, rivers)
# ele.df$RIV_DIST_MIN = apply(x, 1, min)

# ******************************************************************************
#                                       STS
# ******************************************************************************
ele.df <- ele.df %>% nog() %>% 
  relocate(SZN_2, SZN_3, SZN_4, LAT, LON, .after=SEASON) %>% 
  relocate(TOD, .after=DATE.TIME) %>% 
  relocate(CATEGORY, CLASS, COVER, STRUCTURE, VEGETATION, 
           MM.RAIN, MM.RAIN_C, MM.RAIN_7, 
           .after=PROVIDER) %>% 
  dplyr::select(-SEASON)

# STS
setDataPaths('elephant')
save(ele.df, file=procpath("elephant.rdata"))
