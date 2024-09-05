# cleanEle.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanEle.R')
source(here::here('02_scripts', 'utilities.R'))
quickload()
setDataPaths('precipitation')
load(procpath('precipitation.rdata'))
setDataPaths('elephant')
ele.nam <- read.csv(rawpath('nam.eles.fences.csv'))
ele.bots <- read.csv(rawpath('ecoexist.border.fence.pts.csv'))
ele <- rbind(ele.nam, ele.bots)
# ******************************************************************************
#                             Initial looks
# ******************************************************************************
ele$INX <- 1:nrow(ele)
ids <- unique(ele$id)
ele <- ele %>% 
  rename_all(toupper) %>% 
  mutate( 
    # setting up date-time stuff
    DATE.TIME = as.POSIXlt(DATE.TIME, tz="UTC"),
    DTS = DT, 
    DTM = DTS / 60,
    MONTH = month(DATE.TIME),
    YEAR = year(DATE.TIME),
    
    # renaming IDs
    ANIMAL_ID = ID,
    SEX = toupper(SEX),
    ID = match(ANIMAL_ID, ids),
    
    # new burst after days with 0 counts, also new per ID
    DATE = date(DATE.TIME),
    DIFF = DATE - lag(DATE),
    DIFF = ifelse(is.na(DIFF), 1, DIFF),
    BURST = ID*1000 + (cumsum(!DIFF %in% c(0,1))),
    
    # indexing start.count on bursts
    START.COUNT = BURST != lag(BURST),
    START.COUNT = ifelse(is.na(START.COUNT), TRUE, START.COUNT)
  ) 

## add an end count as well.
END.INX = which(ele$START.COUNT) - 1
END.INX = c(END.INX[-1], nrow(ele))
ele$END.COUNT = FALSE
ele$END.COUNT[END.INX] <- TRUE

# distances etc are tied to next point, not the previous... I don't like that.
tolag <- c("DIST", "DX", "DY", "DT", "DTS", "DTM", "R2N",
           "ABS.ANGLE", "REL.ANGLE", "SPEED")
# ele[,tolag] <- lag(ele[,tolag]) 
# ele[ele$START.COUNT,tolag] <- NA
ele$MPS <- ele$DIST / ele$DTS

# set fixrate for collars manually -- I can't figure a way to do this nicely
fillHours <- function(ids, before=NULL, nhours="5 hours") {
  inx <- ele %>% nog() %>% filter(ID %in% ids)
  if (!is.null(before)) {
    before <- as.POSIXlt(before)
    inx <- inx %>% filter(DATE.TIME < before)
  }
  ele$FIXRATE[ele$INX %in% inx$INX] <<- nhours
}

i=i+1
df <- nog(ele) %>% filter(ID == i) %>%
  mutate(DATE.TIME = as.POSIXct(DATE.TIME))
# ggplot(df %>% filter(DTM < 500)) +
#   geom_bar(aes(x=INX, y=DTM, fill=FIXRATE), stat='identity') +
#   ggtitle(paste(i,df$SEX[1], df$COUNTRY[1])) + big.theme
# View(df)

# BOTS notes
# 52-68 start w 4-hour fixes

ele$FIXRATE = "1 hour"
fillHours(1:4,  '2015-08-13 13:00:00', '5 hours')
fillHours(5:6,  '2016-09-09 07:00:00', '5 hours')
fillHours(7:12, '2016-03-22 08:00:00', '5 hours')
fillHours(13:17,'2017-10-11 17:00:00', '5 hours')
fillHours(c(18:19, 22), '2024-01-101 00:00:00', '4 hours')
fillHours(38:43, "2013-11-01 11:00:00", "4 hours")
# started : 2016-05-14 19:00:00
fillHours(52:62, "2016-05-23 11:00:00", "4 hours")
fillHours(63:68, nhours="4 hours")
fillHours(c(63:64, 66:68), "2014-07-15 18:00:00", "8 hours")

# TODO: set some rows to THROW AWAY when running 1hour analyses (data is too frequent)
# ele.df$THROW <- NA
# ID 20 has many of these
# ID 18 and 22 have many of these (4 hours)



# set SEASON based on first rainfall from precipitation.rdata
ele$ISWET <- TRUE
dt <- 121
inx <- which( yday(ele$DATE) > dt )
inx.rain <- match(ele$YEAR+1, wet.start$RAINYEAR)
iswet <- ele$DATE > wet.start$DATE[inx.rain]
ele$ISWET[inx] <- iswet[inx]
ele$SEASON = ifelse(ele$ISWET, "WET", "DRY")

# save to ele.df
ele.df <- ele


# ******************************************************************************
#                                     STATS
# ******************************************************************************
# makeHist <- function(i) {
#   data <- ele.df %>% 
#     filter(ID == i,
#            # DATE.TIME > as.Date('2010-12-01'),
#            # DATE.TIME < as.Date('2011-01-30')
#            ) %>% 
#     mutate(DATE = date(DATE.TIME))
#   hist.data <- data %>% 
#     group_by(DATE, .drop=FALSE) %>% 
#     summarize(n=n()) %>%
#     mutate(FLAG = ifelse(n > 27, "HIGH", ifelse(n<20, "LOW", "AVG")))
#   title = paste0('Elephant ', data$ID[1], ' (', data$SEX[1], ')')
#   cols <- list(AVG='darkgray', HIGH='#08c952', LOW='#f2055c')[unique(hist.data$FLAG)]
#   ggplot() +
#     geom_bar(data=hist.data,
#              mapping=aes(x=DATE, y=n, fill=FLAG),
#              stat="identity") +
#     ggtitle(title) + theme(axis.title.x=element_blank()) +
#     scale_fill_manual(values=cols) + ylab('number of fixes') +
#     scale_x_date(breaks='1 year', date_labels = "%Y") + theme(text=element_text(size=20))
# }
# makeY <- function(i) {
#   data <- ele.df %>% 
#     group_by(ID, SEX) %>% 
#     filter(ID == i, as.n(DTM) > 45) %>% 
#     mutate(
#       BURST = factor(BURST),
#       MONTH = month(DATE.TIME),
#       YEAR = year(DATE.TIME),
#       MONTHYEAR = paste(MONTH, YEAR),
#       MPH = DIST / (as.n(DTM)/60)) 
#   base = ggplot( data=data, 
#           mapping=aes(x=as.Date(DATE.TIME), group=MONTHYEAR, color=SEASON)) + 
#     theme(axis.title.x=element_blank()) + 
#     scale_color_brewer(palette="Dark2", direction=-1) +
#     scale_x_date(breaks='1 year', date_labels = "%Y") + theme(text=element_text(size=20))
#   base
# }

# ----------------------------------
# might have too little data: 26, 27
# eles with higher than average bits of data: 7, 8, 20, 38
# 38-42 have a lot of days with fewer-than-expected pings...
# ----------------------------------


# ******************************************************************************
#                                     CLEANUP BY STATS
# ******************************************************************************

# some data needs to be ignored (collar obviously fell off or died). Use plots above.
findRM <- function(i, start, end=NULL) {
  data <- ele.df %>% filter(ID == i)
  if (is.null(end)) end = data$INX[nrow(data)]
  start : end
}
findRMDate <- function(id, start, end) {
  data <- ele.df %>% 
    filter(ID %in% id, DATE <= as.Date(end), DATE >= as.Date(start))
  data$INX
}
# ele 5 drops at index 95020
inx.5 <- findRM(5, 95020)
# ele 7 drops at index 164495
inx.7 <- findRM(7, 164495)
# ele 19 drops at index 405397
inx.19 <- findRM(19, 405397)
# ele 38 drops at index 612261
inx.38 <- findRM(38, 612261)
# ele 1 needs dec 01 2015 removed
inx.1_4 <- findRMDate(1:4, "2015-12-01", "2015-12-01")
# ele 13 has a weird dip in pings after january 2019
inx.13.1 <- findRMDate(13, "2019-01-09", "2019-01-30")
inx.13.2 <- findRMDate(13, "2019-02-27", "2019-03-08")
# ele 31-37 is missing data in feb 2011
inx.31_37 <- findRMDate(31:37, "2011-02-12", "2011-02-28")
# ele 40 missing a bunch of data in 2015 jan/feb so we will cut off
inx.40 <- findRMDate(40, "2015-01-13", "2015-02-06")

# remove dud data
inx.rm <- c(inx.1_4, inx.5, inx.7, inx.13.1, inx.13.2, inx.19, 
            inx.31_37, inx.38, inx.40)
ele.df <- ele.df[-c(inx.rm),]


# ******************************************************************************
#                            Transforming to spatial data
# ******************************************************************************

ele.df <- ele.df %>% 
  st_as_sf(coords=c('X', 'Y'), crs=32734, remove=FALSE) %>% 
  st_transform(crs = st_crs(khau))

# add latlon
latlon <- st_coordinates(ele.df)
colnames(latlon) <- c('LON', "LAT")
ele.df <- cbind(ele.df, latlon)


# ******************************************************************************
#                                  Final cleanup 
# ******************************************************************************

# reassign INX and pull cols in order
ele.df$INX <- 1:nrow(ele.df)
ele.df$YEAR = year(ele.df$DATE.TIME)
ele.df <- ele.df %>% 
  dplyr::select(INX, ID, BURST, START.COUNT, END.COUNT,
                ANIMAL_ID, SEX,
                DATE.TIME, DATE, YEAR, SEASON, LON, LAT,
                X, Y, DX, DY, DIST, DTM, R2N,
                ABS.ANGLE, REL.ANGLE, MPS,
                FIXRATE, COUNTRY, PROVIDER)
ele.df$DIST[ele.df$END.COUNT] <- NA
ele.df$DTM[ele.df$END.COUNT] <- NA

# ******************************************************************************
#                                       STS
# ******************************************************************************

setDataPaths('elephant')
save(ele.df, file=procpath('elephant.rdata'))
st_write(ele.df, dsn=procpath('elephants.shp'))


