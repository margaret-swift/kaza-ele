# covariateUsePlots.R
# Created 28 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# Statistical models and plots for landscape covariate use values

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/2_eda/eleStats.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(reshape2, nlme)
quickload()
quickload('elephant')


# ******************************************************************************
#                                Covariate modeling
# ******************************************************************************

stats <- nog(ele.khau) %>% 
  group_by(STATE, SZN_4) %>% 
  summarize(evi=median(EVI, na.rm=TRUE),
            evi_sd = sd(EVI, na.rm=TRUE),
            water=mean(waterDist, na.rm=TRUE),
            water_sd = sd(waterDist, na.rm=TRUE))
ggplot( ele.khau ) + 
  geom_histogram(aes(x=RIV_DIST_MIN, fill=STATE), color='gray') + 
  facet_wrap(~SZN_4+STATE, ncol=3, scales="free_y" ) + 
  plot.theme + scale_fill_brewer(palette="Dark2")


# ******************************************************************************
#                                Covariate boxplots
# ******************************************************************************

boxMe <- function(var, df, varname=var, 
                  szn="SZN_2", states=FALSE, cols=NULL) {
  df$YVAR <- df[,var]
  df$SEASON <- df[,szn]
  p=ggplot(df) + 
    geom_boxplot(aes(x=SEASON, y=YVAR, fill=SEASON), color='black') + 
    theme(text=element_text(size=20)) + 
    ylab(varname)
  if (states) p = p + facet_wrap(~STATE)
  if (!is.null(cols)) p = p+scale_fill_manual(values=cols)
  return(p)
}
d
p1 <- boxMe('waterDist', ele.khau %>% filter(STATE != "correlated walk"),
            'distance from \nnatural water (m)') + 
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())
p2 <- boxMe('RIV_DIST_MIN', ele.khau %>% filter(STATE != "correlated walk"),
            'distance from \nrivers (m)')
p1 / p2 + plot_layout(guides="collect")


# ******************************************************************************
#                                Covariate scatterplots
# ******************************************************************************

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
p1 <-scatterMe('waterDist', ele.khau,'distance from natural water (m)') 
p2 <-scatterMe('RIV_DIST_MIN', ele.khau,'distance from rivers (m)')
p1 / p2 + plot_layout(guides="collect")


# ******************************************************************************
#                          Covariate standard deviation plots
# ******************************************************************************

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


# ******************************************************************************
#                               LME Models of covariates
# ******************************************************************************

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
      print( anova(mod) )
    } else if (type == "nostates") {
      mod <- lme(YVAR~SEASON, random=~1|ID, data=df)
      print(mod)
    } else {
      data = df
      mod <- lme(YVAR~SEASON*STATE, random=~1|ID, data=data)
      print(mod$coefficients$fixed)
      print( anova(mod) )
    }
  }
  
  modEach()
  modEach("correlated walk")
  print(d[[1]]$cohen.d)
  modEach("foraging")
  print(d[[2]]$cohen.d)
  modEach("resting")
  print(d[[3]]$cohen.d)
  
  # no states
  message("NO STATES")
  modEach("nostates")
  d = psych::cohen.d(YVAR~SEASON, data=df)
  print(d$cohen.d)
}

library(nlme)
modMe("RIV_DIST_MIN", ele.khau)
modMe("waterDist", ele.khau)
modMe("EVI", ele.khau)


# ******************************************************************************
#                                SSNM powerpoint plots
# ******************************************************************************

# lme model of river distance by season
data <- ele.khau %>% nog() %>% filter(SEX == "F", !is.na(SZN_2)) 
control = lmeControl(msMaxIter = 100)
mod1 <- lme(RIV_DIST_MIN~SEASON, random=~1|ID, data=data, control=control)
summary(mod1)
print( anova(mod1) )
psych::cohen.d(RIV_DIST_MIN~SEASON, data=data)$cohen.d

# lme model of step size by season
distdata <- data %>% filter(DTM <= 60)
mod2 <- lme(DIST~SZN_2 , random=~1|ID, data=distdata)
summary(mod2)
print( anova(mod2) )
psych::cohen.d(DIST~SZN_2, data=distdata)$cohen.d

# boxplot
cols <- c("gray", "#63D4BA")
boxMe('RIV_DIST_MIN', data, varname="", cols=cols )
boxMe('DIST', distdata %>% filter(DIST < 2000), 
      varname="", cols=cols)
