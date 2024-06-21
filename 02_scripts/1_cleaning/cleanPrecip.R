# cleanPrecip.R
# Created 03 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanPrecip.R')
source(here::here('02_scripts', 'utilities.R'))
setDataPaths('precipitation')
precip <- read.csv(rawpath('precipitation_kaza_2010_2022.csv'))


# ******************************************************************************
#                             Wet season start
# ******************************************************************************

dt <- 121
precip.df <- precip %>% 
  mutate(DATE = as.Date(DOY, origin = paste0(YEAR,"-01-01")), 
         INX = 1:nrow(precip),
         RAINYEAR = factor( ifelse(DOY > dt, YEAR+1, YEAR) ),
         RAINDOY = DOY - dt,
         RAINDOY = ifelse(RAINDOY < 0, 365+RAINDOY, RAINDOY)) %>% 
  dplyr::select(INX, DOY, RAINDOY, DATE, RAINYEAR, MM_RAIN) 

# reset rainfall count at the beginning of the dry season
precip.df$MM_RAIN_C <- ave(precip.df$MM_RAIN, 
                           cumsum(precip.df$DOY==dt), 
                           FUN = cumsum)
wt = 1
wet.start <- precip.df %>% 
  group_by(RAINYEAR) %>% 
  filter(MM_RAIN_C >= wt) %>% 
  summarize(RAINDOY = first(RAINDOY)) %>% 
  left_join(precip.df, by=c("RAINYEAR", "RAINDOY"))


# ******************************************************************************
#                             Wet season end
# ******************************************************************************
dt <- 175
precip.df <- precip %>% 
  mutate(DATE = as.Date(DOY, origin = paste0(YEAR,"-01-01")), 
         INX = 1:nrow(precip),
         RAINYEAR = factor( ifelse(DOY > dt, YEAR+1, YEAR) ),
         RAINDOY = DOY - dt,
         RAINDOY = ifelse(RAINDOY < 0, 365+RAINDOY, RAINDOY)) %>% 
  dplyr::select(INX, DOY, RAINDOY, DATE, RAINYEAR, MM_RAIN) 
head(precip.df)

# reset rainfall count at the beginning of the dry season
precip.df$MM_RAIN_C <- ave(precip.df$MM_RAIN, 
                           cumsum(precip.df$DOY==dt), 
                           FUN = cumsum)
wt = 5
wet.end <- precip.df %>% 
  group_by(RAINYEAR) %>% 
  filter(MM_RAIN_C <= wt) %>% 
  summarize(RAINDOY = first(RAINDOY)) %>% 
  left_join(precip.df, by=c("RAINYEAR", "RAINDOY")) %>% 
  mutate(RAINYEAR=as.n( as.c( RAINYEAR ) )-1, 
         RAINDOY = RAINDOY + 300)
NY = nrow(wet.end)
wet.end <- rbind(wet.end[2:NY,], wet.end[NY,])
wet.end$RAINYEAR[NY] = 2023

# ******************************************************************************
#                             plotting
# ******************************************************************************

# rainfall by day
ggplot(mapping=aes(group=RAINYEAR, color=RAINYEAR)) + 
  geom_bar( data=precip.df, 
            mapping=aes(x=RAINDOY, y=MM_RAIN),
            stat="identity") + 
  geom_vline(data=wet.start, aes(xintercept=RAINDOY)) + 
  geom_vline(data=wet.end,
             aes(xintercept=RAINDOY)) + 
  plot.theme + 
  facet_wrap(~RAINYEAR) + 
  xlab('days since start of dry season') + 
  ylab('cumulative rainfall (mm)') 



# cumulative rainfall
ggplot(mapping=aes(group=RAINYEAR, color=RAINYEAR)) + 
  geom_line(data=precip.df, 
            mapping=aes(x=RAINDOY, y=MM_RAIN_C),
            linewidth=1) + 
  geom_vline(data=wet.start, aes(xintercept=RAINDOY)) + 
  plot.theme + 
  facet_wrap(~RAINYEAR) + 
  xlab('days since start of dry season') + 
  ylab('cumulative rainfall (mm)') 

save(wet.start, precip.df,
     file=procpath('precipitation.rdata'))


