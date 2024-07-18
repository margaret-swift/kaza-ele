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
wt = 5
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

rain.theme = big.theme + 
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())
lvl_3 <- c('WET1', 'WET2', 'DRY')
lvl_4 <- c('T1', 'WET', 'T2', 'DRY')
raindat <- precip.df %>%
  mutate(
       MONTH = factor(month(DATE, label=TRUE)),
       YEAR = factor(year(DATE)),
       FAKEDATE = paste0('1899', gsub('\\d{4}', '', DATE)),
       FAKEDATE = as.Date(gsub("02-29", "02-28", FAKEDATE)),
       WEEK = week(FAKEDATE),
       WEEK=ifelse(WEEK==53, 52, WEEK),
       SZN_3 = factor(
         case_when(
           WEEK %in% 1:17 ~ lvl_3[2],
           WEEK %in% 43:53 ~ lvl_3[1],
           .default = lvl_3[3]
         ), levels=lvl_3),
       SZN_4 = factor(
         case_when(
           WEEK %in% 43:48 ~ lvl_4[1],
           WEEK %in% 12:17 ~ lvl_4[3],
           WEEK %in% 18:42 ~ lvl_4[4],
           .default = lvl_4[2]
         ), levels=lvl_4),
       WEEK = factor(WEEK, levels=c(43:53, 1:42))
     ) %>% 
   filter(YEAR %in% 2011:2022)

# rainfall by day of year
ggplot(precip.df) + 
  geom_bar( mapping=aes(x=RAINDOY, y=MM_RAIN),
            stat="identity") + 
  geom_vline(data=wet.start, aes(xintercept=RAINDOY)) + 
  geom_vline(data=wet.end,
             aes(xintercept=RAINDOY)) + 
  rain.theme + 
  facet_wrap(~RAINYEAR, ncol=2) + 
  xlab('days since start of dry season') + 
  ylab('cumulative rainfall (mm)')

# rainfall by day
ggplot(raindat) + 
  geom_bar( mapping=aes(x=DATE, y=MM_RAIN),
            stat="identity") +
  facet_wrap(~YEAR, scales="free", ncol=2) + 
  xlab('date') + ylab('cumulative rainfall (mm)') +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  rain.theme

# rainfall binned by month
raindat %>%
  group_by(YEAR, MONTH) %>% 
  summarize(MM_RAIN = sum(MM_RAIN, na.rm=TRUE)) %>% 
  ggplot() + 
    geom_bar( mapping=aes(x=MONTH, y=MM_RAIN),
              stat="identity") +
    facet_wrap(~YEAR, ncol=2) + 
    xlab('month') + ylab('cumulative rainfall (mm)') +
    rain.theme

# rainfall binned by month, stacked
raindat %>%
  group_by(YEAR, MONTH) %>% 
  summarize(MM_RAIN = sum(MM_RAIN, na.rm=TRUE)) %>% 
  ggplot() + 
  geom_bar( mapping=aes(x=MONTH, y=MM_RAIN, fill=YEAR),
            stat="identity") +
  xlab('month') + ylab('cumulative rainfall (mm)') +
  rain.theme

# rainfall binned by week, stacked
date_lines_3 <- data.frame(
  x=c(0.5, 10.5, 27.5),
  WEEK=c(2.5, 12.5, 29.5), 
  y=36,
  label=c('Oct 30', 'Jan 01', 'Apr 30')
)
date_lines_4 <- data.frame(
  x=c(0.5, 6.5, 21.5, 27.5),
  WEEK=c(2.5, 8.5, 23.5, 29.5), 
  y=36,
  label=c('Oct 30', "Dec 06", 'Mar 01', 'Apr 30')
)
date_lines = date_lines_4
rain_by_week <- raindat %>%
  rename(SEASON=SZN_4) %>% 
  group_by(YEAR, WEEK, SEASON) %>%
  summarize(MM_RAIN = sum(MM_RAIN, na.rm=TRUE)) %>% 
  group_by(WEEK, SEASON) %>% 
  summarize(mu = mean(MM_RAIN, na.rm=TRUE),
            se = std.error(MM_RAIN, na.rm=TRUE))
ggplot(data=rain_by_week, aes(x=WEEK, xend=WEEK)) + 
  geom_vline(data=date_lines, mapping=aes(xintercept=x), 
             linetype='dashed', color='darkgray')  +
  geom_segment(aes(y=mu+se, yend=mu-se, color=SEASON), linewidth=1) +
  geom_point(aes(y=mu, color=SEASON), size=3) +
  geom_label(data=date_lines, 
             mapping=aes(x=WEEK, y=y, label=label))  +
  scale_y_continuous(breaks=seq(0,40,5)) +
  scale_color_brewer(palette="Set2", name="season") +
  xlab('week of year') + ylab('weekly rainfall (mm)') +
  theme(panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(linewidth=1.2)) +
  big.theme + guides(color="none")

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

save(wet.start, wet.end, precip.df,
     file=procpath('precipitation.rdata'))


