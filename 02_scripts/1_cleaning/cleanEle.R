# cleanEle.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanEle.R')
source(here::here('02_scripts', 'utilities.R'))
setDataPaths('elephant')
ele <- read.csv(here(rawpath, 'nam.eles.fences.csv'))

# ******************************************************************************
#                             Initial looks
# ******************************************************************************
ele$INX <- 1:nrow(ele)
ids <- unique(ele$id)
ele.df <- ele %>% 
  rename_all(toupper) %>% 
  mutate(DATE.TIME = as.POSIXlt(DATE.TIME, tz="UTC"),
         ANIMAL_ID = ID,
         START.COUNT = ANIMAL_ID != lag(ANIMAL_ID),
         START.COUNT = ifelse(is.na(START.COUNT), TRUE, START.COUNT),
         SEX = toupper(SEX),
         ID = match(ANIMAL_ID, ids),
         DTS = DT, 
         DTM = difftime(DATE.TIME, lag(DATE.TIME))
         ) %>% 
  dplyr::select(INX, ID, ANIMAL_ID, SEX, START.COUNT, DATE.TIME,
                X, Y, DX, DY, DIST, DTS, DTM, R2N,
                ABS.ANGLE, REL.ANGLE, SPEED,
                COUNTRY, PROVIDER)

ele.pts <- ele.df %>% 
  st_as_sf(coords=c('X', 'Y'), remove=FALSE, crs=32734) %>% 
  st_transform(crs = st_crs(khau))
ele.lines <- ele.pts %>% 
  group_by(ID, ANIMAL_ID, SEX) %>% 
  summarize() %>%
  st_cast("LINESTRING")

## STS
save(ele, ele.df, ele.pts, ele.lines, 
     file=here(procpath, "ele.rdata"))
