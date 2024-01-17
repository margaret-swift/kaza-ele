# determineSeasonalBreaks
# Created 10 Jan 2024
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
pacman::p_load(zoo, here)
i_am('02_scripts/3_modeling/determineSeasonalBreaks.R')
source(here('02_scripts', 'utilities.R'))
setDataPaths('elephant')
load(procpath('ele.rdata'))
setDataPaths('precipitation')
precip <- read.csv(rawpath('precipitation_kaza_2010_2022.csv'))
load(here(outdir, 'hmmLongF.rdata'))


# ******************************************************************************
#                        Calculating cumulative rainfall
# ******************************************************************************
dt <- 200
precip.df <- precip %>% 
  mutate(INX = 1:n(),
         DATE = as.Date(DOY, origin = paste0(YEAR,"-01-01")),
         MM_RAIN_C3 = lag(rollsumr(MM_RAIN, k=3, fill=0)),
         MM_RAIN_C7 = lag(rollsumr(MM_RAIN, k=7, fill=0)),
         MM_RAIN_C10 = lag(rollsumr(MM_RAIN, k=10, fill=0)),
         MM_RAIN_DIFF_C10 = MM_RAIN_C10 - lag(MM_RAIN_C10),
         RAINYEAR = factor( ifelse(DOY > dt, YEAR+1, YEAR) ),
         RAINDOY = DOY - dt,
         RAINDOY = ifelse(RAINDOY < 0, 365+RAINDOY, RAINDOY)
         ) %>% 
  filter(RAINYEAR %in% 2011:2023) %>% 
  dplyr::select(INX, DATE, DOY, RAINYEAR, RAINDOY,
                MM_RAIN, MM_RAIN_C3, MM_RAIN_C7, 
                MM_RAIN_C10, MM_RAIN_DIFF_C10)
  
# ******************************************************************************
#                           seasonal slices
# ******************************************************************************

matchMe <- function(slice) {
  dates <- precip.df %>% mutate(MATCH=paste(RAINYEAR, RAINDOY))
  matchTo <- function(pattern) {
    dates$DATE[match(paste(slice$RAINYEAR, pattern), dates$MATCH)]
  }
  slice$START.DATE = matchTo(slice$START)
  slice$END.DATE = matchTo(slice$END)
  return(slice)
}

# wide slices - START season when rainfall in one 
# ten-day window is greater than 1mm
t1=1
slice.wide <- precip.df %>% 
  filter(MM_RAIN_C10 > t1) %>% 
  group_by(RAINYEAR) %>% 
  summarize(START=min(RAINDOY),
            END=max(RAINDOY)) %>% 
  matchMe()

# narrow slices - only START season when rainfall
# difference from last window is greater than 3mm
t2=3
slice.narrow <- precip.df %>% 
  filter(abs(MM_RAIN_DIFF_C10) > t2) %>% 
  group_by(RAINYEAR) %>% 
  summarize(START=min(RAINDOY),
            END=max(RAINDOY)) %>% 
  matchMe()


# ******************************************************************************
#                                Assigning seasons
# ******************************************************************************

# seasons assigned by cumulative rainfall over past 10 days and split by two, 
#   four, or six seasons per year

slice = slice.wide
szn_2 = slice %>% 
  slice(rep(1:n(), each = 2)) %>% 
  mutate(INX = 1:n(),
         SEASON = rep(c("DRY", "WET"), times=n()/2),
         
         START.DATE=as.c(START.DATE),
         END.DATE=as.c(END.DATE),
         START.DATE = ifelse(SEASON == "DRY", lag(END.DATE), START.DATE),
         END.DATE = ifelse(SEASON == "DRY", lead(START.DATE), END.DATE),
         # set start date for first dry season to 2010-04-25 (average yday 114)
         START.DATE = ifelse(INX == 1, "2010-04-25", START.DATE),
         # set end date for last wet season to 2023-04-25 (average yday 114)
         END.DATE = ifelse(INX == n(), "2023-04-25", END.DATE),
         START.DATE=as.Date(START.DATE),
         END.DATE=as.Date(END.DATE), 
         NDAYS = END.DATE - START.DATE ) %>% 
  dplyr::select(RAINYEAR, SEASON, START.DATE, END.DATE, NDAYS)

szn_4 <- szn_2 %>% 
  slice(rep(1:n(), each = 2)) %>% 
  mutate(
    PART = rep(1:2, times=n()/2),
    SEASON = paste(SEASON, PART, sep="_"),
    NDAYS = NDAYS / 2,
    END.DATE = ifelse(PART==1, START.DATE + NDAYS, END.DATE),
    START.DATE = ifelse(PART==2, lag(END.DATE), START.DATE),
    END.DATE = as.Date(END.DATE),
    START.DATE = as.Date(START.DATE))

szn_6 <- szn_2 %>% 
  slice(rep(1:n(), each = 3)) %>% 
  mutate(
    PART = rep(1:3, times=n()/3),
    SEASON = paste(SEASON, PART, sep="_"),
    NDAYS = NDAYS / 3,
    END.DATE = ifelse(PART==1, START.DATE + NDAYS, END.DATE),
    START.DATE = ifelse(PART==2, lag(END.DATE), START.DATE),
    END.DATE = ifelse(PART==2, START.DATE + NDAYS, END.DATE),
    START.DATE = ifelse(PART==3, lag(END.DATE), START.DATE),
    END.DATE = as.Date(END.DATE),
    START.DATE = as.Date(START.DATE),
    )

szn_2$INTERVAL = interval(szn_2$START.DATE, szn_2$END.DATE)
szn_4$INTERVAL = interval(szn_4$START.DATE, szn_4$END.DATE)
szn_6$INTERVAL = interval(szn_6$START.DATE, szn_6$END.DATE)

# finally adding seasons
assignSzn <- function(d, szn) szn$SEASON[d %within% szn$INTERVAL]
fixSzn <- function(data) {
  fac <- as.factor(as.character(data))
  inx.rm <- which(grepl('c\\(', fac))
  fac[inx.rm] <- NA
  fac <- as.factor(as.character(fac))
}

# add seasons to activity state data
pdata.f$SZN_2 <- sapply(pdata.f$DATE.TIME, function(d) assignSzn(d, szn_2)) %>% fixSzn()
pdata.f$SZN_4 <- sapply(pdata.f$DATE.TIME, function(d) assignSzn(d, szn_4)) %>% fixSzn()
pdata.f$SZN_6 <- sapply(pdata.f$DATE.TIME, function(d) assignSzn(d, szn_6)) %>% fixSzn()

save(hmm.f, pdata.f, data.f, file=here(outdir, 'hmmLongF.rdata'))

inx.rm <- which( ( (pdata.f$SZN_2=="NA") + 
                     (pdata.f$SZN_4=="NA") + 
                     (pdata.f$SZN_6 == "NA") ) > 0 )
pdata <- pdata.f[-inx.rm,]
pdata$SZN_2 <- as.factor(as.character(pdata$SZN_2))
pdata$SZN_4 <- as.factor(as.character(pdata$SZN_4))
pdata$SZN_6 <- as.factor(as.character(pdata$SZN_6))


# ******************************************************************************
#                                 add river distance
# ******************************************************************************

# setDataPaths('geographic')
# load(procpath('geographic.rdata'))
# x = st_distance(pdata.f, rivers)
# pdata.f$RIV_DIST_MIN = apply(x, 1, min)

# ******************************************************************************
#                                     modeling
# ******************************************************************************

pacman::p_load(lme4, lmerTest, emmeans)

mod2 <- lmer(step ~ SZN_2 + STATE + (1 | ID), data=pdata)
mod4 <- lmer(step ~ SZN_4 + STATE + (1 | ID), data=pdata)
mod6 <- lmer(step ~ SZN_6 + STATE + (1 | ID), data=pdata)

# summaries
summary(mod2)
summary(mod4)
summary(mod6)

# mod 6 is the best
AIC(mod2, mod4, mod6)

# Comparing using Anova() because there are more than 2 levels and not nested
# https://www.bookdown.org/rwnahhas/RMPH/mlr-distinctions.html
# model 6 once again comes out on top
car::Anova(mod2, type=3) # Wald Chisq = 256.17  
car::Anova(mod4, type=3) # Wald Chisq = 590.43  
car::Anova(mod6, type=3) # Wald Chisq = 819.42  


# Tukey HS - only one pairing in 6 split is insignificant!
t2 <- lmerTest::difflsmeans(mod2, conf.level=.95)
t4 <- lmerTest::difflsmeans(mod4, conf.level=.95)
t6 <- lmerTest::difflsmeans(mod6, conf.level=.95)
plot(t2)
plot(t4)
plot(t6)

plotSeasonalModels <- function(s) {
  data <- pdata
  if (s == 2) {
    data$SEASON <- data$SZN_2
    mod <- mod2
  } else if (s == 4) {
    data$SEASON <- data$SZN_4
    mod <- mod4
  } else if (s == 6) {
    data$SEASON <- data$SZN_6
    mod <- mod6
  }
  newdata <- data.frame(matrix(nrow=s*3, ncol=3))
  names(newdata) <- c('SEASON', 'STATE', 'ID')
  newdata$STATE <- rep(unique(data$STATE), each=s)
  newdata$SEASON <- rep(levels(data$SEASON), times=3)
  newdata$ID = unique(data$ID)[1]
  
  if (s == 2) {
    x = predict(mod, newdata %>% rename(SZN_2=SEASON))
  } else if (s == 4) {
    x = predict(mod, newdata %>% rename(SZN_4=SEASON))
  } else if (s == 6) {
    x = predict(mod, newdata %>% rename(SZN_6=SEASON))
  }
  newdata$step=x
  ggplot() + 
    geom_point(data=data %>% sample_n(20000),
               mapping=aes(x=SEASON, y=step, color=STATE),
               position='jitter', alpha=0.4) + 
    geom_line(data=newdata,
               mapping=aes(x=SEASON, y=step, group=STATE), size=1) +
    geom_point(data=newdata,
               mapping=aes(x=SEASON, y=step, color=STATE)) + 
    facet_wrap ( ~ STATE, ncol=1, scales="free_y")
}

plotSeasonalModels(2, ys)
plotSeasonalModels(4, ys)
plotSeasonalModels(6, ys)

ggsave(filename=here(outdir, "seasontests", "tukeyhsd.png"))









# ******************************************************************************
#                                     plotting
# ******************************************************************************

# rainfall by day
ggplot(mapping=aes(group=RAINYEAR, color=RAINYEAR)) +
  geom_bar( data=precip.df,
            mapping=aes(x=RAINDOY, y=MM_RAIN_C3),
            stat="identity") +
  # geom_vline(data=wet.START, aes(xintercept=RAINDOY)) +
  plot.theme +
  facet_wrap(~RAINYEAR) +
  xlab('days since START of dry season') +
  ylab('cumulative rainfall (mm)')

# cumulative rainfall
ggplot(mapping=aes(group=RAINYEAR, color=RAINYEAR)) +
  geom_line(data=precip.df,
            mapping=aes(x=RAINDOY, y=MM_RAIN_C),
            linewidth=1) +
  geom_vline(data=wet.START, aes(xintercept=RAINDOY)) +
  plot.theme +
  facet_wrap(~RAINYEAR) +
  xlab('days since START of dry season') +
  ylab('cumulative rainfall (mm)')


ggplot(precip.df, aes(x=RAINDOY, y=MM_RAIN_C10, group=RAINYEAR, color=RAINYEAR)) +
  geom_bar( stat="identity") +
  geom_vline(data=slice.wide, aes(xintercept=START), color='blue') +
  geom_vline(data=slice.wide, aes(xintercept=END), color='red') +
  facet_wrap(~RAINYEAR) +
  guides(color="none") +
  xlab('days since rain-year cutoff') +
  ylab('10-day cumulative rainfall (mm)') +
  plot.theme

ggplot() +
  geom_rect(data=slice.wide,
            mapping=aes(xmin=START.DATE, xmax=END.DATE, ymin=0, ymax=Inf),
            fill='darkblue', alpha=0.1) +
  geom_rect(data=slice.narrow,
            mapping=aes(xmin=START.DATE, xmax=END.DATE, ymin=0, ymax=Inf),
            fill='darkblue', alpha=0.2) +
  geom_bar( data=precip.df,
            mapping=aes(x=DATE, y=MM_RAIN_C10),
            stat="identity") +
  xlab('date') +
  ylab('10-day rainfall (mm)') +
  # ylim(c(0, 15)) +
  plot.theme +
  scale_x_date(date_breaks="year", date_labels='%Y')

plotSzn <- function(s,ys) {
  colors <- c('#fae2a2', '#f5cb5d', '#f0b005',
              '#a8e3da', '#60d6c4', '#02d9b8')
  if (s == 2) {
    data = szn_2
    colors = colors[c(2,5)]
  }else if (s == 4) {
    data = szn_4
    colors = colors[c(1,3,4,6)]
  }else if (s == 6) {
    data = szn_6
  }
  data = data %>% filter(RAINYEAR%in% ys)
  slices <- slice.wide %>% filter(RAINYEAR %in% ys)
  prec = precip.df %>% filter(RAINYEAR%in% ys)
  p=ggplot() +
    geom_rect(data=data,
              mapping=aes(xmin=START.DATE, xmax=END.DATE, ymin=0, ymax=Inf,
                          fill=SEASON)) +
    geom_bar( data=prec,
              mapping=aes(x=DATE, y=MM_RAIN_C10),
              stat="identity") +
    geom_vline(data=slices, aes(xintercept=START.DATE), color='black') +
    geom_vline(data=slices, aes(xintercept=END.DATE), color='black') +
    xlab('date') +
    ylab('10-day rainfall (mm)') +
    plot.theme +
    scale_x_date(date_breaks="year", date_labels='%Y') +
    scale_fill_manual(values = colors)
  return(p)
}

ys <- 2016:2018
p2 <- plotSzn(2, ys) + ylab('') + xlab(element_blank())
p4 <- plotSzn(4, ys) + xlab(element_blank())
p6 <- plotSzn(6, ys) + ylab('')
p2 / p4 / p6


plotSeasonalModels <- function(s) {
  data <- pdata
  if (s == 2) {
    data$SEASON <- data$SZN_2
    mod <- mod2
  } else if (s == 4) {
    data$SEASON <- data$SZN_4
    mod <- mod4
  } else if (s == 6) {
    data$SEASON <- data$SZN_6
    mod <- mod6
  }
  newdata <- data.frame(matrix(nrow=s*3, ncol=3))
  names(newdata) <- c('SEASON', 'STATE', 'ID')
  newdata$STATE <- rep(unique(data$STATE), each=s)
  newdata$SEASON <- rep(levels(data$SEASON), times=3)
  newdata$ID = unique(data$ID)[1]
  
  if (s == 2) {
    x = predict(mod, newdata %>% rename(SZN_2=SEASON))
  } else if (s == 4) {
    x = predict(mod, newdata %>% rename(SZN_4=SEASON))
  } else if (s == 6) {
    x = predict(mod, newdata %>% rename(SZN_6=SEASON))
  }
  newdata$step=x
  ggplot() + 
    geom_point(data=data %>% sample_n(20000),
               mapping=aes(x=SEASON, y=step, color=STATE),
               position='jitter', alpha=0.4) + 
    geom_line(data=newdata,
               mapping=aes(x=SEASON, y=step, group=STATE), size=1) +
    geom_point(data=newdata,
               mapping=aes(x=SEASON, y=step, color=STATE)) + 
    facet_wrap ( ~ STATE, ncol=1, scales="free_y")
}


ggsave(filename=here(outdir, "seasontests", "seasonbreakdown.png"))




colors <- c('#fae2a2', '#f5cb5d', '#f0b005',
            '#a8e3da', '#60d6c4', '#02d9b8')

start.date <- as.Date("2016-04-17")
end.date <- as.Date("2019-05-04")
data <- pdata %>% 
  filter(DATE.TIME > start.date,
         DATE.TIME < end.date) %>%
  mutate(DATE = date(DATE.TIME)) %>% 
  group_by(DATE, SZN_6) %>% 
  summarize(step=mean(step, na.rm=TRUE))
ggplot(data) + 
  geom_bar(mapping=aes(x=DATE, y=step, fill=SZN_6),
             position="stack", stat="identity") + 
  plot.theme + 
  ylab('average step size (m)') +
  scale_fill_manual(values=colors)
           
