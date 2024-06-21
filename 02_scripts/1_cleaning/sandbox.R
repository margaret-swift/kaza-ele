
## CRAWL WRAP SANDBOX
# https://cran.r-project.org/web/packages/momentuHMM/vignettes/momentuHMM.pdf

meta.url <- "https://datarepository.movebank.org/bitstreams/368b24a0-6046-4005-be0c-5a0fa7a24fe3/download"
meta <- read.csv(url(meta.url))

data.url <- "https://datarepository.movebank.org/bitstreams/77e38af2-b5fc-4c10-9e6d-f36b2adf9ea9/download"
rawData <- read.csv(url(data.url))

# select and rename relevant columns
rawData <- rawData[,c(11,3,4,5,6)]
colnames(rawData) <- c("ID","time","lon","lat","temp")

# only keep first track
rawData <- subset(rawData,ID==unique(ID)[1])

# convert times from factors to POSIX
rawData$time <- as.POSIXct(rawData$time,tz="GMT")

# project to UTM coordinates using package rgdal
library(rgdal)
llcoord <- st_as_sf(rawData[,3:4],
                     coords=c('lon', 'lat'))
st_crs(llcoord) <- "+proj=longlat +datum=WGS84"
utmcoord <- st_transform(llcoord, "+proj=utm +zone=30 ellps=WGS84") %>% 
  st_coordinates()

# add UTM locations to data frame
rawData$x <- utmcoord[,'X']
rawData$y <- utmcoord[,"Y"]

dat <- rawData
dat <- dat[-c(10, 20, 15),]

# fit crawl model
crwOut <- crawlWrap(
  obsData=dat, 
  timeStep="hour",
  theta=c(6.855, -0.007), 
  fixPar=c(NA,NA)
  )
View(crwOut$crwPredict)
