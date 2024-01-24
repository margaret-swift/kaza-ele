# eleStats.R
# Created 04 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# Calculate statistics for elephant movement paths. This includes:
# - step size and speed by sex, season, activity state
# - transition matrices for each sex, season
# - step length and turning angles by sex, season, activity state
# - home range size by sex, season

# TODO: 
# - run HMM again with ~SEASON as formula
 
# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/2_eda/eleStats.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(sp, adehabitatHR, reshape2)
setDataPaths('elephant')
load(procpath('ele.rdata'))
# load(here(outdir, 'hmmLongF.rdata'))
# load(here(outdir, 'hmmLongM.rdata'))

setDataPaths('geographic')
load("pdatatmp.rdata")

# ******************************************************************************
#                           General path statistics
# ******************************************************************************

# sort individuals by # days collared
ele.df %>% nog() %>% 
  group_by(ID) %>% 
  summarize(min=min(DATE),
            max=max(DATE),
            days = max-min) %>% 
  arrange(days) %>% 
  print(n=43)

# sort individuals by fixrate type and give count
ele.df %>% nog() %>% 
  group_by(ID, FIXRATE) %>% 
  summarize(n=n()) %>% 
  arrange(FIXRATE) %>% 
  print(n=120)

# average speeds and distances by sex and season
speeds <- ele.df %>% 
  nog() %>% 
  left_join(pdata.f %>% dplyr::select(INX, STATE), by='INX') %>% 
  left_join(pdata.m %>% dplyr::select(INX, STATE), by="INX") %>% 
  mutate(STATE = dplyr::coalesce(STATE.x, STATE.y)) %>% 
  filter(!is.na(STATE))

( speeds %>% 
    group_by(SEX, STATE, SEASON) %>% 
    summarize(DIST.m  = mean(DIST, na.rm=TRUE),
              # DIST.med = median(DIST, na.rm=TRUE),
              DIST.sd  = sd(DIST, na.rm=TRUE),
              # DIST.max = max(DIST, na.rm=TRUE),
              MPS.m   = mean(MPS, na.rm=TRUE),
              # MPS.med = median(MPS, na.rm=TRUE),
              MPS.sd  = sd(MPS, na.rm=TRUE)
              # MPS.max = max(MPS, na.rm=TRUE),
              # n=n())
    )
)


# ******************************************************************************
#                               Transition matrices
# ******************************************************************************

findTransitionMatrix <- function(hmm) {
  vals <- 1:3
  probs.df <- data.frame(matrix(nrow=1, ncol=9))
  names(probs.df) <- paste(rep(vals, each=3), vals, sep=" -> ")
  mat <- hmm$mle$beta
  for (i in vals) {
    types <- paste(i, vals[vals != i], sep=" -> ")
  
    d2 <- exp(mat[,types[1]])
    d3 <- exp(mat[,types[2]])
    sum <- d2 + d3
  
    p2 <- d2 / ( 1 + sum )
    p3 <- d3 / ( 1 + sum )
    p1 <- 1 - ( p2 + p3 )
  
    probs <- data.frame( p1, p2, p3 )
    names(probs) <- c(paste(i, i, sep=" -> "), types)
    probs.df[,names(probs)] <- probs
  }
  probs.df <- unlist(probs.df)
  tm <- matrix(probs.df, nrow=3, byrow=TRUE)
}

# Transition matrices for males and females. 
# TODO: run HMM for seasons and grab transition matrices for both sexes and 
#   both seasons
tm.f <- findTransitionMatrix(hmm.f)
tm.m <- findTransitionMatrix(hmm.m)



# ******************************************************************************
#                                 Home ranges
# ******************************************************************************

# Set up data to find home ranges for each i-s-y
ref.table <- ele.df %>% nog() %>%  
  mutate(YEAR = year(DATE.TIME)) %>% 
  group_by(ID, SEX, SEASON, YEAR) %>% 
  summarize(NPOINTS=n(), 
            DATE.START = min(DATE.TIME),
            DATE.END = max(DATE.TIME),
            NDAYS = (date(DATE.END) - date(DATE.START)),
            NDAYS = as.numeric(NDAYS)) %>% 
  filter(NDAYS > 100) %>% #has to have enough data
  mutate(id = paste(ID, YEAR, SEASON, sep="_"), 
         ID = as.c(ID))

# Pull area of each individual-season-year
sp <- ele.df %>% nog() %>% 
  mutate(
    YEAR = year(DATE.TIME),
    id = paste(ID, YEAR, SEASON, sep="_"))
sp <- SpatialPointsDataFrame(coords=sp[,c('X', 'Y')], data=sp)[,c('id')]
# kud <- kernelUD(sp, h='href')
# hr <- lapply(kud[1:10], function(e)
  # getverticeshr(e, percent=95, unin='m', unout='m2'))
mcp.df <- mcp(sp, percent=95) %>% as.data.frame()

ref.table %>% 
  left_join(mcp.df, by='id') %>% 
  group_by(SEX, SEASON) %>% 
  summarize(AREA_M = mean(area, na.rm=TRUE),
            AREA_SD = sd(area, na.rm=TRUE))





# ******************************************************************************
#                           MATCHING PATTERNS TO STATISTICS
# ******************************************************************************

# Elephant movements trace along fences and channel along omiramba



# Pinch points across roads and at rivers	



# Habitual/repeated movement to and from water sources, especially artificial waterholes		



# Deflection/permeability differences by sex and boundary type (fence [and fence type], road, river)



# Elephants move away from water sources in the wet season (expansion contraction)



# 80% of roan crossings were during the wet season. Does this happen for elephant too?	



# Elephant attraction to some areas with higher quality resources?	