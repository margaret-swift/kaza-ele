# m4.R
# Created 16 Oct 2023
# Margaret Swift <mes114@duke.edu>

# Cullen, J. A., Poli, C. L., Fletcher, R. J., & Valle, D. (2022).
#   Identifying latent behavioural states in animal movement with M4,
#   a nonparametric Bayesian method. Methods in Ecology and Evolution,
#   13, 432–446. https://doi.org/10.1111/2041-210X.13745

# GitHub page: https://joshcullen.github.io/bayesmove/
# Vignettes: https://joshcullen.github.io/bayesmove/articles/index.html

# NOTES
# - okay so this is working pretty well until the piece at the end of
# this doc when we try to establish breakpoints. Consistently only like 6-12
# are being establishd, which is not nearly enough. How do I figure this out?

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

pacman::p_load(here, bayesmove, furrr, dplyr)
i_am('02_scripts/3_modeling/m4.R')
source(here('02_scripts', 'utilities.R'))
makeDataPaths('elephant')
load(here(procpath, 'ele.rdata'))

# ******************************************************************************
#                                 PREPPING DATA
# ******************************************************************************
ggplot(ele.df %>% mutate(BURST = factor(BURST)) %>% filter(ID <= 4)) +
  geom_line(aes(x=INX, y=DX, color=BURST, group=BURST)) +
  facet_wrap(~ID, nrow=4, scales="free_x")

good <- c(1001, 2004, 3006, 4009)
# good <- c(1001)
data <- ele.df %>%
  nog() %>%
  dplyr::filter(BURST %in% good,
                FIXRATE == "1 hour") %>%
  dplyr::select(ID, BURST, X, Y, DATE.TIME) %>%
  rename(date = DATE.TIME, id = BURST) %>%
  as.data.frame()
# data <- data[1:5000,]

tracks <- prep_data(dat = data,
                    coord.names = c("X", "Y"),
                    id = "id")
tracks <-
  round_track_time(
    dat = tracks,
    id = "id",
    int = 3600,
    tol = 180,
    time.zone = "UTC",
    units = "secs"
  )
inx.na <- is.na(tracks$angle)
tracks <- tracks[!inx.na, ]

# Create list from data frame
tracks.list <- df_to_list(dat = tracks, ind = "id")

# Filter observations
tracks_filt.list <- filter_time(dat.list = tracks.list, int = 3600)

# View sample of results
# head(tracks_filt.list[[3]])

# Check that only observations at 1 hour time intervals are retained per ID
purrr::map(tracks_filt.list, ~ n_distinct(.$dt))

# Define bin number and limits for turning angles
angle.bin.lims = seq(from = -pi, to = pi, by = pi / 4)  #8 bins
angle.bin.lims[c(1, 9)] <- c(-3.3, 3.3)

# Define bin number and limits for step lengths
dist.bin.lims = quantile(tracks[tracks$dt == 3600, ]$step,
                         # c( seq(0, 1, length.out=5)), na.rm=T)  #5 bins
                         c(0, 0.25, 0.50, 0.75, 0.90, 1), na.rm = T)  #5 bins

# Assign bins to observations
tracks_disc.list <- map(
  tracks_filt.list,
  discrete_move_var,
  lims = list(dist.bin.lims, angle.bin.lims),
  varIn = c("step", "angle"),
  varOut = c("SL", "TA")
)

# # Histograms of step length and turning angle
# hist(tracks$step, col='lightblue', breaks=200,
#      xlim=c(0, 8000))
# abline(v=dist.bin.lims, col='purple', lty=2, lwd=2)
# hist(tracks$angle, col='pink', breaks=50)
# abline(v=angle.bin.lims, col='purple', lty=2, lwd=2)



# ******************************************************************************
#                                 SEGMENTING DATA
# ******************************************************************************

# Only retain id and discretized step length (SL) and turning angle (TA) columns
tracks.list2 <- map(tracks_disc.list,
                    subset,
                    select = c(id, SL, TA))

set.seed(1)

# Define hyperparameter for prior distribution
alpha <- 1

# Set number of iterations for the Gibbs sampler
ngibbs <- 1e5

# Set the number of bins used to discretize each data stream
nbins <- c(length(dist.bin.lims),
           length(angle.bin.lims)) - 1

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan(multisession, workers = 3)  #run all MCMC chains in parallel
#refer to future::plan() for more details
dat.res <- segment_behavior(
  data = tracks.list2,
  ngibbs = ngibbs,
  nbins = nbins,
  alpha = alpha
)

future::plan(future::sequential)  #return to single core


# Trace-plots for the number of breakpoints per ID
# traceplot(data = dat.res, type = "nbrks")

# Trace-plots for the log marginal likelihood (LML) per ID
# traceplot(data = dat.res, type = "LML")

# Determine MAP for selecting breakpoints
NR = nrow(dat.res$LML)
MAP.est = vector(mode = "numeric", length = NR)
for (i in 1:NR) {
  x = get_MAP(dat = dat.res$LML[i, ], nburn = 1000)
  MAP.est[i] <- x
  message(i, ' : ', x)
}

# select breakpoints
brkpts <- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)
brkpts <- brkpts[2:ncol(brkpts)]

# How many breakpoints estimated per ID?
apply(brkpts, 1, function(x)
  length(purrr::discard(x, is.na)))

# assign segments to data
mytracks = rbi::flatten(tracks_disc.list)
mytracks$tseg <- NA
mytracks$INX <- data.table::rowid(mytracks$id)
myrows <- match(mytracks$id, rownames(brkpts))
i = 57926
for (i in 1:nrow(mytracks)) {
  j <- myrows[i]
  brks <- purrr::discard(brkpts[j, ], is.na)
  mytracks$tseg[i] <- sum(!(mytracks$INX[i] <= brks)) + 1
}
tracks.seg <-
  mytracks # assign_tseg(dat = tracks_disc.list, brkpts = brkpts)

# # plotting breakpoints
# my.plot_breakpoints <- function(i) {
#   data <- tracks_disc.list[[i]]
#   data$INX <- 1:nrow(data)
#   breaks <- melt(purrr::discard(brkpts[i,][-1], is.na))
#   (
#     p1 <- ggplot(data=data, mapping=aes(x=INX, y=step)) +
#       geom_line(color='orange') +
#       geom_vline(data=breaks,
#                  mapping=aes(xintercept = value),
#                  linewidth=1)
#   )
#   (
#     p2 <- ggplot(data=data, mapping=aes(x=INX, y=angle)) +
#       geom_line(color='darkgreen') +
#       geom_vline(data=breaks,
#                  mapping=aes(xintercept = value),
#                  linewidth=1)
#     )
#
#   p1 + p2 + plot_layout(nrow=2)
# }
# my.plot_breakpoints(3)



# ******************************************************************************
#                                 CLUSTER SEGMENTS
# ******************************************************************************

tracks.seg2 <- tracks.seg[, c("id", "tseg", "SL", "TA")]

# Summarize observations by track segment
obs <- summarize_tsegs(dat = tracks.seg2, nbins = nbins)

# Prepare for Gibbs sampler
ngibbs <- 1e4  #number of MCMC iterations for Gibbs sampler
nburn <- ngibbs / 2  #number of iterations for burn-in
nmaxclust <-
  max(nbins) - 1  #one fewer than max number of bins used for data streams
ndata.types <- length(nbins)  #number of data types

# Priors
gamma1 <- 0.1
alpha <- 0.1

set.seed(1000)
# Run LDA model
res <- cluster_segments(
  dat = obs,
  gamma1 = gamma1,
  alpha = alpha,
  ngibbs = ngibbs,
  nmaxclust = nmaxclust,
  nburn = nburn,
  ndata.types = ndata.types
)

# Check traceplot of log likelihood
# plot(res$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")

# Extract proportions of behaviors per track segment
theta.estim <- extract_prop(
  res = res,
  ngibbs = ngibbs,
  nburn = nburn,
  nmaxclust = nmaxclust
)

# Convert to data frame for ggplot2
theta.estim_df <- theta.estim %>%
  as.data.frame() %>%
  pivot_longer(
    .,
    cols = 1:all_of(nmaxclust),
    names_to = "behavior",
    values_to = "prop"
  ) %>%
  modify_at("behavior", factor)
levels(theta.estim_df$behavior) <- 1:nmaxclust

# Plot results
# ggplot(theta.estim_df, aes(behavior, prop)) +
#   geom_boxplot(fill = "grey35",
#                alpha = 0.5,
#                outlier.shape = NA) +
#   geom_jitter(color = "grey35",
#               position = position_jitter(0.1),
#               alpha = 0.3) +
#   labs(x = "\nBehavior", y = "Proportion of Total Behavior\n")

# Extract bin estimates from phi matrix
behav.res <-
  get_behav_hist(
    dat = res,
    nburn = nburn,
    ngibbs = ngibbs,
    nmaxclust = nmaxclust,
    var.names = c("Step Length", "Turning Angle")
  )

# Plot histograms of proportion data
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  scale_fill_manual(
    values = c(
      "#21908CFF",
      "#440154FF",
      "#FDE725FF",
      "grey35",
      "grey35",
      "grey35",
      "grey35"
    ),
    guide = "none"
  ) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")


# Reformat proportion estimates for all track segments
behav.names.4 = c("Encamped", "Transit", "Foraging Short", "Foraging Long")
behav.names.3 = c("Encamped", "Transit", "Foraging")
behav.order.4 = c(1,3,4,2)
behav.order.3 = c(1,3,2)

theta.estim.long <- expand_behavior(
  dat = tracks.seg,
  theta.estim = theta.estim,
  obs = obs,
  nbehav = 3,
  behav.names = behav.names.3,
  behav.order = behav.order.3
  # behav.names = c("ARS", "Encamped", "Transit"),
  # behav.order = c(2, 1, 3)
)
unique(theta.estim.long$behavior)

# Plot results
# ggplot(theta.estim.long) +
#   geom_area(
#     aes(x = date, y = prop, fill = behavior),
#     color = "black",
#     linewidth = 0.25,
#     position = "fill"
#   ) +
#   labs(x = "\nTime", y = "Proportion of Behavior\n") +
#   scale_fill_viridis_d("Behavior") +
#   facet_wrap( ~ id)

# Convert segmented dataset into list
tracks.list <- df_to_list(dat = tracks.seg, ind = "id")

# Merge results with original data
tracks.out <- assign_behavior(
  dat.orig = mytracks,
  dat.seg.list = tracks.list,
  theta.estim.long = theta.estim.long,
  behav.names = behav.names.3
)

# Map dominant behavior for all IDs
data <-tracks.out[tracks.out$id == '4009',]
ggplot(data=tracks.out, mapping=aes(x=x, y=y)) +
  geom_path(  color = "gray60", linewidth = 0.25) +
  geom_point(
    mapping=aes(fill = behav),
    pch=21,
    alpha = 0.5
  ) +
  scale_fill_viridis_d("Behavior", na.value = "grey50") +
  facet_wrap( ~ id, scales = "free", ncol = 2)


# Proportion encamped
plotTrack <- function(id, type) {
  data <- tracks.out[tracks.out$id == id,]
  firstlast <- data[c(1,nrow(data)),c('x', 'y', 'INX')]
  names(data)[names(data) == type] <- 'TARGET'
  keytitle <- paste0("Proportion\n", type)
  p <- ggplot() +
    geom_path( data=data, aes(x, y, color = TARGET),
               linewidth = 0.5, alpha = 0.7) +
    geom_point(data = firstlast, aes(x, y, shape=factor(INX)), size=3) +
    scale_color_distiller(keytitle, palette = "Spectral", na.value = "grey50") +
    labs(x = "Easting", y = "Northing", title = type) +
    facet_wrap( ~ id, scales = "free")
  p
}
plotAll <- function(id) {
  data <- tracks.out[tracks.out$id == id,]
  data$behav = as.c(data$behav)
  data$behav[is.na(data$behav)] = data$behav[which(is.na(data$behav))-1]
  data$behav[is.na(data$behav)] = "UNKNOWN"
  data$behav <- factor(data$behav)
  data$flag = data$behav == lag(data$behav)
  data$DATE = date(data$date)
  
  firstlast <- data[c(1,nrow(data)),c('x', 'y', 'INX')]
  
  p <- ggplot() +
    geom_path( data=data, 
               mapping=aes(x, y, color=behav, group=interaction(behav, DATE)),
               linewidth = 0.5, alpha = 0.7) +
    geom_point(data = firstlast, aes(x, y, shape=factor(INX)), size=3) +
    scale_color_manual(keytitle, values=c('red', 'blue', 'purple', 'gray50')) +
    labs(x = "Easting", y = "Northing", title = type) +
    facet_wrap( ~ id, scales = "free")
  p
}
plotAll('1001')
