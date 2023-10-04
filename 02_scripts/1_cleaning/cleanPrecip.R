# cleanPrecip.R
# Created 03 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanPrecip.R')
source(here::here('02_scripts', 'utilities.R'))
setDataPaths('precipitation')
precip <- read.csv(here(rawpath, 'precipitation_kaza_2010_2022.csv'))
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
head(precip.df)

wt = 5
wet.start <- precip.df %>% 
  group_by(RAINYEAR) %>% 
  filter(MM_RAIN_C >= wt) %>% 
  summarize(RAINDOY = first(RAINDOY)) %>% 
  left_join(precip.df, by=c("RAINYEAR", "RAINDOY"))
  
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
     file=here(procpath, 'precipitation.rdata'))


