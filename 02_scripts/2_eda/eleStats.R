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
load(procpath('eleKhau.rdata'))
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
#                                Covariate modeling
# ******************************************************************************


stats <- nog(ele.khau) %>% group_by(STATE, SZN_4) %>% 
  summarize(evi=median(EVI, na.rm=TRUE),
            evi_sd = sd(EVI, na.rm=TRUE),
            water=mean(waterDist, na.rm=TRUE),
            water_sd = sd(waterDist, na.rm=TRUE)
  )

ggplot( ele.khau ) + 
  geom_histogram(aes(x=RIV_DIST_MIN, fill=STATE), color='gray') + 
  facet_wrap(~SZN_4+STATE, ncol=3, scales="free_y" ) + 
  plot.theme + scale_fill_brewer(palette="Dark2")


boxMe <- function(var, df, varname=var) {
  df$YVAR <- df[,var]
  p=ggplot(df) + 
    geom_boxplot(aes(x=SZN_4, y=YVAR, fill=SZN_4), color='black') + 
    plot.theme + 
    # scale_fill_manual(values=c('#D5A684', '#9CADFF')) +
    # scale_fill_brewer(palette="Dark2") + 
    facet_wrap(~STATE) +
    ylab(varname)
  return(p)
}

boxMe("EVI", ele.khau)
p1 <- boxMe('waterDist', ele.khau %>% filter(STATE != "correlated walk"),
            'distance from \nnatural water (m)') + 
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())
p2 <- boxMe('RIV_DIST_MIN', ele.khau %>% filter(STATE != "correlated walk"),
            'distance from \nrivers (m)')
p1 / p2 + plot_layout(guides="collect")


scatterMe <- function(var, df, varname=var) {
  df$YVAR <- df[,var]
  year(df$DATE) <- 2010 
  s <- df %>% 
    filter(!is.na(SZN_6)) %>% 
    group_by(DATE) %>% 
    summarize(dist = mean(YVAR), sd=sd(YVAR))
  p=ggplot(s) + 
    geom_point(aes(x=DATE, y=dist), color='black', size=2) + 
    geom_segment(aes(x=DATE, xend=DATE,
                     y=dist-sd, yend=dist+sd), color='gray') +
    scale_x_date(date_breaks="1 month", date_labels="%m") +
    ylab(varname) + xlab('date (months)') +
    plot.theme
  p
}

scatterMe("EVI", ele.khau)
p1 <-scatterMe('waterDist', ele.khau,'distance from natural water (m)') 
p2 <-scatterMe('RIV_DIST_MIN', ele.khau,'distance from rivers (m)')
p1 / p2 + plot_layout(guides="collect")

sdMe <- function(var, df, varname=var) {
  df$YVAR <- df[,var]
  year(df$DATE) <- 2010 
  s <- df %>% 
    filter(!is.na(SZN_6)) %>% 
    group_by(DATE) %>% 
    summarize(dist = mean(YVAR), sd=sd(YVAR))
  p=ggplot(s) + 
    geom_line(aes(x=DATE, y=sd), color='black', linewidth=1) + 
    scale_x_date(date_breaks="1 month", date_labels="%m") +
    ylab(varname) + xlab('date (months)') +
    plot.theme
  p
}

p0 <- sdMe('DIST', ele.khau, 's.d. in step size')
p1 <-sdMe('waterDist', ele.khau,'s.d. in dist from \nnatural water (m)') 
p2 <-sdMe('RIV_DIST_MIN', ele.khau,'s.d. in dist from \nrivers (m)')
p0 / p1 / p2 + plot_layout(guides="collect")

# running models
mySummary <- function(mod) {
  summary(mod)
}
modMe <- function(var, df) {
  df$YVAR <- df[,var]
  d = psych::cohen.d(YVAR~SEASON*STATE, data=df)
  
  modEach <- function(type="all") {
    message("SUMMARY FOR ", toupper(type))
    if (type != "all") {
      data = df %>% filter(STATE==type)
      mod <- lme(YVAR~SEASON, random=~1|ID, data=data)
      c <- round(mod$coefficients$fixed[['SEASONWET']], 2)
      message('coefficient: ', c)
    } else {
      data = df
      mod <- lme(YVAR~SEASON*STATE, random=~1|ID, data=data)
      print(mod$coefficients$fixed)
    }
    print( anova(mod) )
  }
  
  modEach()
  modEach("correlated walk")
  print(d[[1]]$cohen.d)
  modEach("foraging")
  print(d[[2]]$cohen.d)
  modEach("resting")
  print(d[[3]]$cohen.d)
}

modMe("waterDist", ele.khau)
modMe("RIV_DIST_MIN", ele.khau)
modMe("EVI", ele.khau)


## Habitat choice model
mod <- lme(YVAR~SEASON*STATE, random=~1|ID, data=data)


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