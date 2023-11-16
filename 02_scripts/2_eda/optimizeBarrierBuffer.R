# determineBarrierBehavior.R
# Created 09 November 2023
# Margaret Swift <mes114@duke.edu>

# Xu, W., Dejid, N., Herrmann, V., Sawyer, H. & Middleton, A. D. 
#   Barrier Behaviour Analysis (BaBA) reveals extensive effects of fencing on 
#   wide‐ranging ungulates. J Appl Ecol 58, 690–698 (2021). 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2664.13806
#   https://github.com/wx-ecology/BaBA/

# Path straightness is the ratio between the displacement distance and the 
# accumulated step length of a trajectory, ranging from 0 (sinuous) to 1 (straight) 

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

# devtools::install_github("wx-ecology/BaBA", force=TRUE)
# rpath <- "C:/Users/mes473/OneDrive - Cornell University/Documents/R_Packages/BaBA/R"
# source(file.path(rpath, "BaBA_fix.R"))
pacman::p_load(sp, BaBA, here)
source(here("02_scripts", "utilities.R"))
setDataPaths('elephant')
load(procpath('ele.rdata'))

# ******************************************************************************
#                                   FUNCTIONS
# ******************************************************************************

# By applying fence buffer distances every 10 meters from 50m - 150m, we
# define the optimal buffer as the distance at which the number of quick cross
# events begins to level off (< 1 % increase).
prepareQC <- function(points, fences, d.list) {
  qc.df <- data.frame(d=d.list, nqc=0, ntot=0, perc=0)
  
  for (i in 1:nrow(qc.df)) {
    d <- d.list[i]
    message("BaBA for ", d, "m")
    tic()
    baba.out <- BaBA(
      animal = points,
      barrier = fences,
      d = d,
      b_time = 10,
      p_time = 100,
      w = 10,
      interval = 2,
      units = "hours",
      tolerance = tolerance,
      max_cross = 3
    )
    toc()
    classes <- baba.out$classification 
    if (nrow(classes)) {
      qc.df$nqc[i] <- sum(grepl('Cross', classes$eventTYPE))
      qc.df$ntot[i] <- nrow(classes)
    }
  }
  qc.df$perc <- qc.df$nqc / qc.df$ntot
  return(qc.df)
}

# ******************************************************************************
#                               OPTIMIZE BUFFER d 
# ******************************************************************************

# Data points mutation
points <- ele.df %>% 
  mutate(Animal.ID = ID, date=DATE.TIME) %>% 
  filter(FIXRATE == '1 hour') %>% 
  filter(SEX == "M") %>% 
  mutate(date = as.POSIXct(date))

# for some reason, st_union needs us to project to UTM, NAD83 seems to work??
points_proj <- st_transform(points, 3578)
fences_proj <- st_transform(fences, 3578)

# run on projected data
qc.df.1h3 <- prepareQC(points_proj, fences_proj, seq(2000, 2500, by=100))
qc.df.all <- rbind(qc.df.all, qc.df.1h3)

# ******************************************************************************
#                                   PLOTTING
# ******************************************************************************

pt.1h <- data.frame(x=1000, y=qc.df.all$perc[qc.df.all$d==1000] )
ggplot() +
  geom_point(data=qc.df.all, aes(x=d, y=perc), color='black') +
  geom_line( data=qc.df.all, aes(x=d, y=perc), color='black') +
  geom_point(data=pt.1h, mapping=aes(x=x, y=y), color="red", size=3) +
  xlab("distance (m)") +
  ylab("proportion 'cross'") + plot.theme 
ggsave(here('03_output', 'eda', 'optimize_d.png'), 
       width=10, height=3.5)
