# trajectoryStats.R
# Created 03 Sept 2024
# Margaret Swift <margaret.swift@cornell.edu>

# Part of simulation validation portion of project
# Vignette from trajr: https://cran.rstudio.com/web/packages/trajr/vignettes/trajr-vignette.html

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
# 
# packages and HERE
pacman::p_load(here, trajr)
i_am('02_scripts/5_validation/trajectoryStats.R')

# utilities
source(here('02_scripts', 'utilities.R'))

# loading data
quickload('elephant')
dat <- ele.df %>% 
  rename_all(tolower) %>% 
  mutate(times=date.time)
setOutPath('elephants')

# ******************************************************************************
#                             FIRST PLOTS
# ******************************************************************************

# pull first burst
x <- dat %>% 
  filter(burst == 1000) %>% 
  select('x', 'y', 'times')
trj <- TrajFromCoords(x[1:1000,])
plot(trj, main="BURST - 1000")
points(trj, draw.start.pt = FALSE, pch = 16, col = "black", cex = 0.7)

# smoothed trajectory
smoothed <- TrajSmoothSG(trj, p = 3, n = 31)
lines(smoothed, col = "#FF0000A0", lwd = 2)
legend("topright", c("Original track", "Smoothed track"), 
       lwd = c(1, 2), lty = c(1, 1), col = c("black", "red"), inset = 0.01)


# ******************************************************************************
#                             ANALYZE SPEED
# ******************************************************************************

# The functions TrajVelocity and TrajAcceleration estimate velocity and 
#  acceleration as vectors at each point along a trajectory. Velocity is 
#  calculated using first-order finite differences, and acceleration uses 
#  second-order finite differences. The TrajDerivatives function calculates
#  speed and change in speed along a Trajectory. If the trajectory is noisy, 
#  it should be smoothed before the derivatives are calculated (see Smoothing 
#  trajectories, but take note of the caveats described there).

# Calculate speed and acceleration
derivs <- TrajDerivatives(smoothed)

# Plot change-in-speed and speed
plot(derivs$acceleration ~ derivs$accelerationTimes, type = 'l', col = 'red', 
     yaxt = 'n',
     xlab = 'Time (s)',
     ylab = expression(paste('Change in speed (', m/s^2, ')')))
axis(side = 2, col = "red")
lines(derivs$speed ~ derivs$speedTimes, col = 'blue')
axis(side = 4, col = "blue")
mtext('Speed (m/s)', side = 4, line = 3)
abline(h = 0, col = 'lightGrey')

# Calculate slow intervals
intervals <- TrajSpeedIntervals(trj, slowerThan = 2000)
print(intervals)
plot(intervals)


# ******************************************************************************
#                                 ANALYZE SINUOSITY
# ******************************************************************************

# Various methods to measure the straightness, or conversely, tortuosity, of 
#  trajectories are available within trajr. The simplest is 𝐷/𝐿, where 𝐷 is 
#  the distance from the start to the end of the trajectory, and 𝐿 is the 
#  length of the trajectory. This straightness index is calculated by the 
#  function TrajStraightness, and is a number ranging from 0 to 1, where 1 
#  indicates a straight line. Benhamou (2004) considers the straightness index 
#  to be a reliable measure of the efficiency of a directed walk, but 
#  inapplicable to random trajectories. Batschelet (1981) considers this 
#  straightness index to be an approximation of r, which is the length of 
#  the mean vector of turning angles after rediscretizing to a constant 
#  step length. r can be calculated by calling Mod(TrajMeanVectorOfTurningAngles(trj)), 
#  assuming trj is a Trajectory with constant step length.
# 
# Within this section, most of the figures show the appropriate statistic plotted 
#  for multiple randomly generated correlated random walks which vary in the 
#  standard deviations of angular errors, (higher values of which produce more 
#  tortuous trajectories).
#  

# Calculate straightness (D/L)and Sinuosity for trajectory at a range of timesteps
createStraightnessStepsDF <- function(trj, NS=100) {
  steps <- seq(0, max(trj$displacementTime)/2, length.out=NS)[2:NS]
  df <- data.frame(time=steps, strt=NA, sinus=NA, emax=NA)
  for (i in 1:(NS-1)) {
    dt <- steps[i]
    y <- trj %>% filter(displacementTime < dt)
    strt <- TrajStraightness(y)
    sinus <- TrajSinuosity2(y)
    emax <- TrajEmax(y)
    
    if (length(strt) && !is.nan(strt)) df$strt[i] <- strt
    if (length(sinus) && !is.nan(sinus)) df$sinus[i] <- sinus
    if (length(emax) && !is.nan(emax) && emax>0) df$emax[i] <- emax
  }
  df
}

characteriseTrajectory <- function(trj) {
  
  # Step and turn
  sl <- TrajStepLengths(trj)
  mean_sl <- mean(sl, na.rm=TRUE)
  sd_sl <- sd(sl, na.rm=TRUE)
  max_sl <- max(sl, na.rm=TRUE)
  ta <- TrajAngles(trj)
  mean_ta <- mean(ta, na.rm=TRUE)
  sd_ta <- sd(ta, na.rm=TRUE)
  
  # Measures of speed
  derivs <- TrajDerivatives(trj)
  mean_speed <- mean(derivs$speed)
  sd_speed <- sd(derivs$speed)
  
  # Measures of straightness
  sinuosity <- TrajSinuosity2(trj)
  Emax <- TrajEmax(trj)
  straightness <- TrajStraightness(trj)
  
  # Return a list with all of the statistics for this trajectory
  data.frame(mean_sl = mean_sl,
             max_sl = max_sl,
             sd_sl = sd_sl,
             mean_ta = mean_ta,
             sd_ta = sd_ta,
             mean_speed = mean_speed,
             sd_speed = sd_speed,
             sinuosity = sinuosity,
             Emax = Emax,
             straightness=straightness
  )
}


# calculate straightness over multiple bursts
NS <- 100
ids <- unique(dat$id)
bursts <- unique(dat$burst[dat$id %in% ids])
strt.df <- NULL
stat.df <- NULL
for (b in bursts) {
  message('running for ele burst ', b, '.')
  x <- dat %>% filter(burst == b)
  
  # function to add info to df's that is the same across time steps
  addInfo <- function(df) {
    m <- x[1,]
    df %>% mutate(id=m$id, fixrate=m$fixrate, sex=m$sex, country=m$country,
                  burst=b) %>% 
      relocate(id, burst, sex, country, fixrate)
  }
 
  # only run if burst is large enough
  if (nrow(x) > 2000) {
    # create trajectory
    trj <- x %>% select('x', 'y', 'times') %>% TrajFromCoords()
    trj <- trj[2:(nrow(trj)-1),]
    message('traj created')
    
    # calculate straightness and add to main DF
    df <- createStraightnessStepsDF(trj, NS) %>% addInfo()
    if (is.null(strt.df)) {
      strt.df <- df
    } else {
      strt.df <- rbind(strt.df, df)
    }
    
    # characterize trajectories
    df <- characteriseTrajectory(trj) %>% addInfo()
    if (is.null(stat.df)) {
      stat.df <- df
    } else {
      stat.df <- rbind(stat.df, df)
    }
  } else { message('burst too short! moving on.') }
}

# make everything factors
fixFactors <- function(df, cols=c('id', 'burst', 'sex')) mutate_at(df, .vars=cols, .funs=as.factor)
strt.df <- fixFactors(strt.df)
stat.df <- fixFactors(stat.df)

# get means
base.df <- strt.df %>% 
  filter(fixrate == "1 hour")
  
bounds = base.df %>% 
  mutate(time = round(time/5)*5) %>% 
  group_by(sex, country, fixrate, time) %>% 
  summarize(mu_sin=mean(sinus, na.rm=TRUE),
            sd_sin=sd(sinus, na.rm=TRUE),
            mu_st=mean(strt, na.rm=TRUE),
            sd_st=sd(strt, na.rm=TRUE)) %>% 
  mutate(min_sin=mu_sin-sd_sin, max_sin=mu_sin+sd_sin,
         min_st=mu_st-sd_st, max_st=mu_st+sd_st)

rcolor <- 'yellow'
ralpha <- 0.2
balpha <- 0.5
ccolor <- c('black','red'); scolor <- c('purple', 'black')
base <- ggplot(base.df, aes(x=time, group=burst, color=sex)) + 
  theme_bw() + big.theme + 
  facet_wrap(~sex+country, scales='free_y', nrow=2) +
  scale_x_continuous(limits=c(0, 50)) + 
  scale_color_manual(values=scolor)
  # scale_color_manual(values=c('black', 'purple'))
p1 <- base +
  geom_ribbon(data=bounds, 
              aes(ymin=min_sin, ymax=max_sin, group=sex), 
              fill=rcolor, alpha=ralpha, color=rcolor) +
  geom_line(data=bounds, aes(y=mu_sin, group=sex),
             color=rcolor, linewidth=1) +
  geom_line(aes(y=sinus), alpha=balpha, linewidth=1) + 
  ylab('sinuosity')
p2 <- base +  
  geom_ribbon(data=bounds, 
              aes(ymin=min_st, ymax=max_st, group=sex), 
              fill=rcolor, alpha=ralpha, color=rcolor) +
  geom_line(data=bounds, aes(y=mu_st, group=sex),
             color=rcolor, linewidth=1) +
  geom_line(aes(y=strt), alpha=balpha, linewidth=1) + 
  ylab('straightness')
p1
fname <- outpath('elephant_movement_statistics', 'sinuosity_by_sex_cntry.png')
ggsave(file=fname, width=15, height=8)

p2
fname <- outpath('elephant_movement_statistics', 'straightness_by_sex_cntry.png')
ggsave(file=fname, width=15, height=8)

# normalize stats
cols <- c('mean_sl', 'max_sl', 'mean_speed', 'sinuosity')
normalize <- function(x, na.rm = T) (x / max(x, na.rm = T))
stat.norm <- stat.df %>% 
  mutate_at(cols, normalize) %>% 
  mutate(group = paste(country, sex)) %>% 
  select(id, burst, sex, country, all_of(cols)) %>% 
  melt()

ggplot(stat.norm, aes(x=country, y=value, fill=sex)) + 
  geom_boxplot(position="dodge") +
  facet_wrap(~variable, nrow=2) +
  theme_bw() + big.theme + 
  scale_fill_manual(values=scolor)
#
fname <- outpath('elephant_movement_statistics', 'means_by_sex_cntry.png')
ggsave(file=fname, width=7, height=7)











