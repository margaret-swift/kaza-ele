# findRelativeSpaceUse.R
# Created 3 May 2024
# Margaret Swift <margaret.swift@cornell.edu>

# amt vignette: 
#   https://cran.r-project.org/web/packages/amt/vignettes/p4_SSF.html
# interpreting SSF results: 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2656.13441


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/habitat_analyses/runHabitatSelection.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(amt, reshape2)

quickload() #loads common spatial features

# load step selection dat
setOutPath('habitat_selection')
load(outpath('fitmodelBN.rdata'))


# ******************************************************************************
#                         COMPUTE RELATIVE SPACE USE
# ******************************************************************************

# find relative use of one habitat type over another
findUse <- function(hsf, v) {
  cs <- coefficients(hsf) %>% t() %>% as.data.frame()
  use <- exp(sum(cs * c(1, v)))
  use
}
findUseAvail <- function(hsf, v, val) {
  avail <- table(ssf.df$cover_class)[[val]] / nrow(ssf.df)
  use <- findUse(hsf, v)
  return(use * avail)
}
vals <- c(mean(ssf.df$evi, na.rm=TRUE))


# bare, crop, open, sparse, closed
val.bare   <- c(vals, c(0,0,0,0))
val.crop   <- c(vals, c(1,0,0,0))
val.open   <- c(vals, c(0,1,0,0))
val.sparse <- c(vals, c(0,0,1,0))
val.closed <- c(vals, c(0,0,0,1))

use.bare <- findUse(hsf, val.bare)
use.crop <- findUse(hsf, val.crop)
use.open <- findUse(hsf, val.open)
use.sparse <- findUse(hsf, val.sparse)
use.closed <- findUse(hsf, val.closed)

ua.bare <- findUseAvail(hsf, val.bare, "bare/water")
ua.crop <- findUseAvail(hsf, val.crop, 'cropland')
ua.open <- findUseAvail(hsf, val.open, 'open')
ua.sparse <- findUseAvail(hsf, val.sparse, 'sparse')
ua.closed <- findUseAvail(hsf, val.closed, 'closed')

# make it a matrix!
top.u <- t(c(use.bare, use.crop, use.open, use.sparse, use.closed))
top.a <- t(c(ua.bare, ua.crop, ua.open, ua.sparse, ua.closed))


makeMat <- function(top) {
  left = t(1/top)
  mat <- left %*% top
  row.names(mat) <- colnames(mat) <- levels(data$cover_class)
  mat
}

heatPlot <- function(vals, inx) {
  # make matrix
  mat <- makeMat(vals)
  m <- t(round(mat[inx,inx], 2))
  diag(m) <- NA
  
  # set up data for plotting
  data <- expand.grid(X=rownames(m), Y=colnames(m))
  data$ratio <- NA
  for (i in 1:nrow(data)) data$ratio[i] <- m[data$X[i], data$Y[i]]
  data$X <- factor(data$X, levels=rev(colnames(m)))
  data$logratio <- log(data$ratio)
  ggplot(data=data, aes(x=X, y=Y)) + 
    geom_tile(aes(fill=logratio), alpha=0.8) + 
    geom_text(aes(label=ratio), size=10) + 
    geom_text(aes(label=paste(X,Y, sep=":"),
                  alpha=X==Y), size=5, vjust=-3) + 
    geom_text(aes(label=ratio), size=10) + 
    scale_fill_distiller(palette="PiYG", type="div", 
                         na.value="lightgray", name="log(ratio)") + 
    theme(text=element_text(size=24),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_blank()) +
    guides(alpha="none") +
    scale_x_discrete(position = "top") +
    scale_alpha_manual(values=c(1,0.35)) 
}
heatPlot(top.u, 3:5) + ggtitle('landcover use ratios, all else equal')
heatPlot(top.a, 3:5) + ggtitle('landcover use ratios, given availability')


# ******************************************************************************
#                              RUN SIMPLE SSF
# ******************************************************************************

m0 <- ssf.df |> fit_clogit(case_ ~ cover + strata(step_id_))
m1 <- ssf.df |> fit_clogit(case_ ~ hillslope + evi + strata(step_id_))
m2 <- ssf.df |> fit_clogit(case_ ~ hillslope + evi + waterdist + cover_class + strata(step_id_))
summary(m0)
summary(m1)
summary(m2)

tbl <- lands.meta[,c('category', 'cover')]
lr.reclass = terra::classify(lr, tbl)
plot(lr.reclass, 
     col=c("gray", "darkgray", "lightgreen", "lightblue", "#af7e7f"))




# ******************************************************************************
#                              HSF FOR ALL LANDCOVER TYPES
# ******************************************************************************

# habitat selection function
hsf <- glm(case_ ~ category * state, 
           data = ssf.df, weight = w, family = binomial)

# summary
summary(hsf)




# ******************************************************************************
#                                     PLOT STATS
# ******************************************************************************

statsMe <- function(var, season=NULL) {
  df = as.data.frame(ssf.df)
  df$VAR = df[,var]
  
  if (is.null(season)) {
    stats = df %>% 
      dplyr::select(case_, VAR) %>% 
      group_by(case_) %>%
      table()
  } else {
    df$SEASON = df[,season]
    a <- df %>% 
      dplyr::select(case_, VAR, SEASON) %>% 
      group_by(case_, SEASON) %>%
      table() %>% 
      reshape2::melt() %>% 
      dcast(formula=case_+SEASON~VAR)
    n <- a[,1:2] %>% mutate(case_=ifelse(case_, "used", "available"))
    stats <- a %>% dplyr::select(-case_,-SEASON) %>% as.matrix()
  }
  
  NR=nrow(stats); NC=ncol(stats);
  tot <- rowSums(stats)
  num = matrix(stats, nrow=NR)
  den = matrix(rep(tot, times=NC), nrow=NR)
  x = data.frame(t(num / den))
  
  # fix col names
  if (!is.null(season)) colnames(x) <- paste(n$case_, n$SEASON, sep="_")
  else colnames(x) = c('used', 'available')
  x$class = factor(colnames(stats), levels=colnames(stats))
  
  return( x )
}
plotStats <- function(x=NULL, xlab=NULL, season=NULL) {
  if (is.null(x)) x = statsMe(var, season)
  data <- reshape2::melt(x, variable.name='case', value.name='freq')
  by_season = length(names(x)) > 3
  
  if (by_season)  {
    v <- data$case
    data$case <- gsub("_.*", "", v)
    data$season <- gsub(".*_", "", v)
  } 
  
  if (is.null(xlab)) xlab <- "class"
  colors <- c('#D01C8B', '#4DAC26')
  colors <- c( 'black', 'darkgray')
  p = ggplot(data) +
    geom_bar(mapping=aes(x=class, y=freq, fill=case), 
             stat="identity", position="dodge") + 
    xlab(xlab) + ylab('frequency of use') +
    plot.theme + 
    scale_fill_manual(values=colors) 
  if (by_season) p <- p + facet_wrap(~season, ncol=1)
  p
}

# cover type, veg, water distance bin
var <- 'cover_class'
season <- "szn_2"
cover.tab <- statsMe(var)

var <- 'veg_class'
veg.tab <- statsMe(var, season)

# var <- 'waterbin'
# water.tab <- statsMe(var, season)

var <- 'category'
season <- "szn2_"
class.tab <- statsMe(var)

# plot used vs available coefficients sensu Fieberg et al 2021 Fig. 1
# https://onlinelibrary.wiley.com/doi/pdf/10.1111/1365-2656.13441
cover.plot <- plotStats(cover.tab, "cover type")
veg.plot <- plotStats(veg.tab, "vegetative class")
water.plot <- plotStats(water.tab, "distance from surface water")
class.plot <- plotStats(class.tab, "landcover class")

set1 <- c(5,1,3)
set2 <- c(5,2,4)
cover.tab[,set1]
cover.tab[,set2]
veg.tab[,set1]
veg.tab[,set2]
water.tab[,set1]
water.tab[,set2]

cover.plot
veg.plot
water.plot
class.plot