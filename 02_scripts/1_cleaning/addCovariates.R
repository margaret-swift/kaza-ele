# addCovariates.R
# Created 04 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/addCovariates.R')
source(here::here('02_scripts', 'utilities.R'))
pacman::p_load(suncalc, zoo, lwgeom)

setDataPaths('geographic')
load(procpath('geographic.rdata'))
# load(procpath('landcover_rasters.rdata'))

setDataPaths('precipitation')
load(procpath('precipitation.rdata'))
load(procpath('season_slices.rdata'))

setDataPaths('elephant')
load(procpath('ele.rdata'))
# 
# # ******************************************************************************
# #                               Precipitation
# # ******************************************************************************
# 
# # add precip in the last week
# precip.df$MM.RAIN_7 <- rollsum(precip.df$MM_RAIN, k=7, fill=NA, align="right")
# 
# p.inx <- match(ele.df$DATE, precip.df$DATE)
# ele.df$MM.RAIN <- precip.df$MM_RAIN[p.inx]
# ele.df$MM.RAIN_C <- precip.df$MM_RAIN_C[p.inx]
# ele.df$MM.RAIN_7 <- precip.df$MM.RAIN_7[p.inx]
# 
# 
# # ******************************************************************************
# #                               Landcover type
# # ******************************************************************************
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
# 
# # ******************************************************************************
# #                             Adding time of day
# # ******************************************************************************
# 
# data <- ele.df %>% nog() %>% 
#   rename(lat=LAT, lon=LON, date.time=DATE.TIME) %>% 
#   mutate(date = as.Date(as.POSIXct(date.time)))
# sun.times <- getSunlightTimes(data=data[,c('date', 'lat', 'lon')])
# findDiff <- function(e) difftime( e, data$date.time )
# 
# dawn <- difftime (sun.times$dawn, data$date.time ) %>% as.n()
# sunriseEnd <- difftime (sun.times$sunriseEnd, data$date.time ) %>% as.n()
# sunset <- difftime (sun.times$sunset, data$date.time ) %>% as.n()
# night <- difftime (sun.times$night, data$date.time ) %>% as.n()
# df <- data.frame(dawn, sunriseEnd, sunset, night)
# 
# # account for before dawn
# signs <- sign(df) == -1
# sums <- rowSums(signs)
# sums[sums == 0] <- 4
# time.names <- c('DAWN', 'DAY', 'DUSK', "NIGHT")
# TOD <- time.names[sums]
# ele.df$TOD <- TOD


# ******************************************************************************
#                             Adding distance from water
# ******************************************************************************
ele.df$SZN_2 <- sapply(ele.df$DATE.TIME, function(d) assignSzn(d, szn_2)) %>% fixSzn()
ele.df$SZN_4 <- sapply(ele.df$DATE.TIME, function(d) assignSzn(d, szn_4)) %>% fixSzn()
ele.df$SZN_6 <- sapply(ele.df$DATE.TIME, function(d) assignSzn(d, szn_6)) %>% fixSzn()


# ******************************************************************************
#                             Adding distance from water
# ******************************************************************************
x = st_distance(ele.df, rivers)
ele.df$RIV_DIST_MIN = apply(x, 1, min)

# ******************************************************************************
#                                       STS
# ******************************************************************************
ele.df <- ele.df %>% 
  relocate(SZN_2, SZN_4, SZN_6, LAT, LON, .after=SEASON) %>% 
  relocate(TOD, .after=DATE.TIME) %>% 
  relocate(MM.RAIN, MM.RAIN_C, MM.RAIN_7, 
           LC_CLASS, LC_CATEG, RIV_DIST_MIN, .before=geometry)

# STS
setDataPaths('elephant')
save(ele.df, ele, file=procpath("ele.rdata"))
