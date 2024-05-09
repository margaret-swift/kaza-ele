# runStepSelection.R
# Created 13 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# amt vignette: 
#   https://cran.r-project.org/web/packages/amt/vignettes/p4_SSF.html
# interpreting SSF results: 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2656.13441


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/habitat_analyses/runIntegratedStepSelection.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(amt, reshape2)
setDataPaths('geographic')
load(procpath('geographic.rdata'))
load(procpath('stepSelectionParamsISSF.rdata'))


# ******************************************************************************
#                              RUN SIMPLE iSSF
# ******************************************************************************
ssf_dat <- ssf.df %>% mutate(cos_ta_=cos(ta_), log_sl_=log(sl_))

m1 <- ssf_dat %>% 
  fit_issf(case_ ~ cover_end + waterdist_end + szn_2 + 
             sl_ + log_sl_ + cos_ta_ + 
             szn_2 : sl_ + szn_2:log_sl_ +
             strata(inx), model = TRUE)
summary(m1)



# ******************************************************************************
#                              CALCULATE RSS
# ******************************************************************************

msl_ = mean(ssf_dat$sl_, na.rm=TRUE)
mwd_ = mean(ssf_dat$waterdist_end, na.rm=TRUE)

# data.frame for s1
s1.wet <- data.frame(
  cover_end = factor("sparse", levels = levels(ssf_dat$cover_end)),
  sl_ = msl_,
  waterdist_end=seq(0, 10000, 500),
  szn_2 = factor("WET", levels=levels(ssf_dat$szn_2)),
  log_sl_ = log(msl_),
  cos_ta_ = 1,
  inx=606918 )

# data.frame for s2
s2.wet <- data.frame(
  cover_end = factor("sparse", levels = levels(ssf_dat$cover_end)),
  waterdist_end=mwd_,
  szn_2 = factor("WET", levels=levels(ssf_dat$szn_2)),
  sl_ = msl_,
  log_sl_ = log(msl_),
  cos_ta_ = 1,
  inx=606918 )

s1.dry <- data.frame(
  cover_end = factor("sparse", levels = levels(ssf_dat$cover_end)),
  sl_ = msl_,
  waterdist_end=seq(0, 10000, 500),
  szn_2 = factor("DRY", levels=levels(ssf_dat$szn_2)),
  log_sl_ = log(msl_),
  cos_ta_ = 1,
  inx=606918 )

# data.frame for s2
s2.dry <- data.frame(
  cover_end = factor("sparse", levels = levels(ssf_dat$cover_end)),
  waterdist_end=mwd_,
  szn_2 = factor("DRY", levels=levels(ssf_dat$szn_2)),
  sl_ = msl_,
  log_sl_ = log(msl_),
  cos_ta_ = 1,
  inx=606918 )

lr1 <- log_rss(m1, x1=s1, x2=s2)
lr1$df


# ******************************************************************************
#                              CALCULATE RSS
# ******************************************************************************


# Plot using ggplot2
ggplot(lr1$df, aes(x = waterdist_end_x1, y = exp(log_rss))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = mwd_, linetype = "dashed", color = "red", alpha=0.5) +
  theme_bw() + plot.theme + 
  xlab('distance from water') + 
  ylab('log-RSS vs Mean waterdist')

lr2_ci_wet <- log_rss(m1, s1.wet, s2.wet, ci = "se", 
                       ci_level = 0.95)

lr2_ci_dry<- log_rss(m1, s1.dry, s2.dry, ci = "se", 
                       ci_level = 0.95)

ggplot() +
  geom_ribbon(data=lr2_ci_wet$df,
              mapping=aes(x = waterdist_end_x1, y = log_rss, ymin = lwr, ymax = upr), 
              linetype = "dashed", 
              color = "black", fill = "gray80", alpha = 0.5) +
  geom_ribbon(data=lr2_ci_dry$df,
              mapping=aes(x = waterdist_end_x1, y = log_rss, ymin = lwr, ymax = upr), 
              linetype = "dashed", 
              color = "black", fill = "gray80", alpha = 0.5) +
  geom_line(data=lr2_ci_wet$df, 
            mapping=aes(x = waterdist_end_x1, y = log_rss),
            color='blue',
            size = 1) +
  geom_line(data=lr2_ci_dry$df,
            mapping=aes(x = waterdist_end_x1, y = log_rss),
            color='brown',
            size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = mwd_, linetype = "dashed", color = "red", alpha=0.5) +
  theme_bw() + plot.theme + 
  xlab('distance from water') + 
  ylab('log-RSS vs Mean waterdist')

ggplot(lr2_ci_dry$df, aes(x = waterdist_end_x1, y = log_rss)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              linetype = "dashed", 
              color = "black", fill = "gray80", alpha = 0.5) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = mwd_, linetype = "dashed", color = "red", alpha=0.5) +
  theme_bw() + plot.theme + 
  xlab('distance from water') + 
  ylab('log-RSS vs Mean waterdist')



# ******************************************************************************
#                             PLOTTING COVER TYPE PARAMS
# ******************************************************************************

updateGamma <- function(mod, covertype) {
  slvar <- paste0('sl_:cover_start', covertype)
  logslvar <- paste0('log_', slvar)
  coefs <- mod$model$coefficients
  beta_sl = coefs["sl_"];  beta_log_sl = coefs["log_sl_"];
  if (slvar %in% names(coefs)) {
    beta_sl = beta_sl + coefs[slvar]
    beta_log_sl = beta_log_sl + coefs[logslvar]
  }
  sl <- update_gamma(
    dist = mod$sl_,
    beta_sl = beta_sl,
    beta_log_sl = beta_log_sl)
  sl
}
updateVM <- function(mod, covertype) {
  costavar <- paste0('cos_ta_:cover_start', covertype)
  coefs <- mod$model$coefficients
  beta_cos_ta = coefs["cos_ta_"]
  if (slvar %in% names(coefs)) beta_cos_ta = beta_cos_ta + coefs[costavar]
  ta <- update_vonmises(
    dist = m3$ta_, 
    beta_cos_ta = beta_cos_ta)
  ta
}

# set up model
m3 <- ssf_dat %>% 
  fit_issf(case_ ~ cover_end + waterdist_end + 
             sl_ + log_sl_ + cos_ta_ + 
             cover_start:(sl_ + log_sl_ + cos_ta_) +
             strata(inx), model = TRUE)
summary(m3)

# iSSF step-length distribution updates
sparse_sl <- updateGamma(m3, "sparse")
open_sl <- updateGamma(m3, "open")
closed_sl <- updateGamma(m3, "closed")

# iSSF turning-angle distribution updates
sparse_ta <- updateVM(m3, "sparse")
open_ta <- updateVM(m3, "open")
closed_ta <- updateVM(m3, "closed")

# data.frame for plotting
plot_sl <- data.frame(x = rep(NA, 100))

# x-axis is sequence of possible step lengths
plot_sl$x <- seq(from = 0, to = 400, length.out = 100)

# y-axis is the probability density under the given gamma distribution
# Forest
plot_sl$open <- dgamma(x = plot_sl$x, 
                         shape = open_sl$params$shape,
                         scale = open_sl$params$scale)
# Grass
plot_sl$sparse <- dgamma(x = plot_sl$x, 
                        shape = sparse_sl$params$shape,
                        scale = sparse_sl$params$scale)
# Wet
plot_sl$closed <- dgamma(x = plot_sl$x, 
                      shape = closed_sl$params$shape,
                      scale = closed_sl$params$scale)

# Pivot from wide to long data
plot_sl <- plot_sl %>% 
  pivot_longer(cols = -x)

# Plot
p1<-ggplot(plot_sl, aes(x = x, y = value, color = factor(name))) +
  geom_line(size = 1) +
  scale_color_manual(name = "Land-use",
                     breaks = c("closed", "open", "sparse"),
                     values = c("forestgreen", "wheat", "blue")) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()


#data.frame for plotting
plot_ta <- data.frame(x = rep(NA, 100))

# x-axis is sequence of possible step lengths
plot_ta$x <- seq(from = -pi, to = pi, length.out = 100)

# y-axis is the probability density under the given gamma distribution
# closed
plot_ta$closed <- circular::dvonmises(x = plot_ta$x, 
                                      kappa = closed_ta$params$kappa,
                                      mu = 0)
# open
plot_ta$open <- circular::dvonmises(x = plot_ta$x, 
                                     kappa = open_ta$params$kappa,
                                     mu = 0)
# sparse
plot_ta$sparse <- circular::dvonmises(x = plot_ta$x, 
                                   kappa = sparse_ta$params$kappa,
                                   mu = 0)

# Pivot from wide to long data
plot_ta <- plot_ta %>% 
  pivot_longer(cols = -x)

# Plot
p2 <- ggplot(plot_ta, aes(x = x, y = value, color = factor(name))) +
  geom_line(size = 2) +
  scale_color_manual(name = "Land-use",
                     breaks = c("closed", "sparse", "open"),
                     values = c("lightgreen", "#af7e7f", "lightblue")) +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = c(expression(-pi, -pi/2, 0, pi/2, pi))) +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  theme_bw() + theme(text=element_text(size=24))
p2

combined <- p1 + p2 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
